"""Connectome audit: can a deletion-based network map from this data be trusted?

Written for anyone running node-deletion analyses on structural connectomes,
which includes connectomic surgical planning. It answers three questions about a
specific dataset rather than about the method in general:

  1. How much of the deletion damage is just node strength? If a measure's
     ranking is degree, it can be reproduced by summing a row of the matrix, and
     nothing is gained by computing it.

  2. Does the ranking survive a change of edge-weighting convention? A map that
     reorders when the same tractography is weighted differently cannot support
     a clinical decision, however good its p-values look.

  3. If it does not survive, what fixes it? Compressing the weight distribution
     turns out to help more than thresholding, and keeps every edge.

The calibration behind the prediction, measured on HCP streamline connectomes
across a grid of tail exponents and densities (run_mechanism.py):

    tail ratio (max/median of nonzero weights)     expected map stability
        ~8700                                          +0.32
         ~900                                          +0.46
         ~190                                          +0.64
          ~90                                          +0.63
          ~30                                          +0.69
          ~14                                          +0.77
          ~10                                          +0.84

Spearman between log tail ratio and stability was -0.917 across that grid, and
density added essentially nothing once the tail was known (partial +0.267).

    python3 audit.py --input connectomes.npy
    python3 audit.py --input data/averageConnectivity_Fpt.csv --log10
"""

import argparse
import os

import numpy as np
from scipy.stats import spearmanr

from controllability import average_controllability, spectral_scale
from lesion import global_efficiency
from run_confounds import residualize
from run_weightings import WEIGHTINGS, reweight

# Tail ratio thresholds, from the calibration grid above.
TAIL_SAFE = 30.0
TAIL_RISKY = 200.0


def _square_numeric(arr):
    """Is this a numeric array with at least one square face?"""
    if arr is None or not hasattr(arr, "dtype"):
        return False
    if arr.dtype.kind not in "fiu":  # float, signed int, unsigned int
        return False
    shape = tuple(arr.shape)
    return len(shape) >= 2 and any(
        shape[i] == shape[j] and shape[i] > 1
        for i in range(len(shape)) for j in range(i + 1, len(shape)))


def _parcel_count(shape):
    """Length of the repeated (square) axis, i.e. the parcel count."""
    for i in range(len(shape)):
        for j in range(i + 1, len(shape)):
            if shape[i] == shape[j] and shape[i] > 1:
                return shape[i]
    return -1


def _better(cand_shape, best_shape, max_nodes):
    """Prefer the largest parcellation that stays within `max_nodes`.

    Multi-scale releases ship the same connectomes at several resolutions, so
    taking the biggest array picks a 1000-parcel scale and makes the audit take
    hours. Anything at or under the cap beats anything over it, and above the
    cap the smallest is the least bad.
    """
    c, b = _parcel_count(cand_shape), _parcel_count(best_shape)
    if b < 0:
        return True
    c_ok, b_ok = c <= max_nodes, b <= max_nodes
    if c_ok != b_ok:
        return c_ok
    return c > b if c_ok else c < b


def _from_mat(path, max_nodes=400):
    """Pull the connectome out of a .mat, both v7.3 and older.

    A .mat usually carries labels and ids alongside the matrices, so the
    candidate has to be numeric and square rather than merely the first array
    encountered. Nested object arrays, which older MATLAB structs produce, are
    walked recursively.
    """
    try:
        import h5py
        with h5py.File(path, "r") as f:
            best, best_shape = None, ()
            for k in f:
                if k.startswith("#"):
                    continue
                d = f[k]
                if not hasattr(d, "shape") or not _square_numeric(d):
                    continue
                if _better(tuple(d.shape), best_shape, max_nodes):
                    best, best_shape = np.array(d), tuple(d.shape)
            if best is None:
                raise ValueError(f"no square numeric matrix in {path}")
            return best
    except OSError:
        pass

    import scipy.io as sio
    m = sio.loadmat(path)
    best, best_shape = None, ()

    def walk(v):
        nonlocal best, best_shape
        v = np.asarray(v)
        # Structured dtypes must be checked before object: a dtype like
        # [('SC', 'O'), ('FC', 'O')] is not equal to object, so testing for
        # object first silently skips every MATLAB struct.
        if v.dtype.names:
            for name in v.dtype.names:
                for item in np.ravel(v[name]):
                    walk(item)
            return
        if v.dtype == object:
            for item in np.ravel(v):
                walk(item)
            return
        if _square_numeric(v) and _better(tuple(v.shape), best_shape, max_nodes):
            best, best_shape = v, tuple(v.shape)

    for k, v in m.items():
        if not k.startswith("__"):
            walk(v)
    if best is None:
        raise ValueError(f"no square numeric matrix in {path}")
    return best


def load_any(path, log10=False, max_nodes=400):
    """Accept .npy, .csv/.txt, or MATLAB .mat, returning (subjects, N, N)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npy":
        A = np.load(path)
    elif ext in (".csv", ".txt", ".tsv"):
        A = np.genfromtxt(path, delimiter="," if ext == ".csv" else None)
    elif ext == ".mat":
        A = _from_mat(path, max_nodes=max_nodes)
    else:
        raise ValueError(f"unsupported file type {ext}")

    A = np.asarray(A, dtype=float)
    if log10:
        A = np.power(10.0, A)
    if A.ndim == 2:
        A = A[None, ...]
    elif A.ndim == 3 and A.shape[0] not in A.shape[1:]:
        pass  # already (subjects, N, N)
    elif A.ndim == 3:
        # (N, N, subjects) is the common MATLAB layout
        if A.shape[0] == A.shape[1]:
            A = np.transpose(A, (2, 0, 1))
    return A


def clean(M):
    M = np.array(M, dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    return (M + M.T) / 2.0


def tail_ratio(M):
    off = M[~np.eye(M.shape[0], dtype=bool)]
    nz = off[off > 0]
    return float(nz.max() / np.median(nz)) if nz.size else float("nan")


def power_transform(M, alpha):
    """Raise nonzero weights to a power, compressing the tail scale-free.

    Unlike log10(1 + w) this behaves the same whether weights are streamline
    counts in the millions or normalized densities below 1, which matters
    because log10(1 + w) is nearly linear for small w and compresses nothing.
    """
    out = M.copy()
    nz = out > 0
    out[nz] = out[nz] ** alpha
    return out


def recommend_alpha(mats, target_tail, grid=None):
    """Smallest compression that brings the tail under `target_tail`.

    Returns the largest alpha (least distortion) meeting the target, or None.
    """
    grid = grid if grid is not None else np.arange(0.95, 0.04, -0.05)
    for alpha in grid:
        tails = [tail_ratio(power_transform(M, alpha)) for M in mats]
        if np.median(tails) < target_tail:
            return float(alpha)
    return None


def damage_sweep(A, measure="ac", max_parcels=None):
    """Per-parcel deletion damage. `max_parcels` subsamples for speed."""
    n = A.shape[0]
    targets = np.arange(n)
    if max_parcels and n > max_parcels:
        targets = np.linspace(0, n - 1, max_parcels).astype(int)

    if measure == "ac":
        ref = spectral_scale(A)
        intact = average_controllability(A, ref_scale=ref)
        out = np.zeros(len(targets))
        for j, i in enumerate(targets):
            keep = np.setdiff1d(np.arange(n), [i])
            out[j] = intact[keep].mean() - average_controllability(
                A[np.ix_(keep, keep)], ref_scale=ref).mean()
    else:
        base = global_efficiency(A)
        out = np.array([base - global_efficiency(A, [i]) for i in targets])
    return out, A.sum(axis=1)[targets]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="connectome file")
    ap.add_argument("--log10", action="store_true",
                    help="stored values are log10; exponentiate first")
    ap.add_argument("--subjects", type=int, default=8,
                    help="how many subjects to test (default 8)")
    ap.add_argument("--max-nodes", type=int, default=400,
                    help="when a file holds several parcellation scales, prefer"
                         " the largest at or under this size")
    ap.add_argument("--max-parcels", type=int, default=200,
                    help="subsample parcels above this count, for speed")
    args = ap.parse_args()

    A = load_any(args.input, log10=args.log10, max_nodes=args.max_nodes)
    n_subj = min(args.subjects, A.shape[0])
    n_parcel = A.shape[1]
    mats = [clean(A[i]) for i in range(n_subj)]

    print(f"\nCONNECTOME AUDIT: {args.input}")
    print(f"  {A.shape[0]} matrices of {n_parcel} parcels, auditing {n_subj}\n")

    # 1. Shape of the data
    dens = np.mean([np.mean(M[~np.eye(n_parcel, dtype=bool)] > 0) for M in mats])
    tails = np.array([tail_ratio(M) for M in mats])
    print("1. WEIGHT DISTRIBUTION")
    print(f"   density                     {dens:.3f}")
    print(f"   tail ratio (max/median)     {np.median(tails):.0f}")
    if np.median(tails) < TAIL_SAFE:
        regime = "compressed. Convention sensitivity should be low."
    elif np.median(tails) < TAIL_RISKY:
        regime = "moderate. Some convention sensitivity expected."
    else:
        regime = "HEAVY. Expect strong convention sensitivity."
    print(f"   regime                      {regime}\n")

    # 2. Degree redundancy
    print("2. IS THE DAMAGE MEASURE JUST NODE DEGREE?")
    for measure, label in (("ac", "average controllability"),
                           ("ge", "global efficiency")):
        rs = []
        for M in mats:
            d, st = damage_sweep(M, measure, args.max_parcels)
            rs.append(spearmanr(d, st)[0])
        r = float(np.mean(rs))
        verdict = ("redundant with degree" if r > 0.85
                   else "largely degree" if r > 0.6 else "carries extra signal")
        print(f"   {label:<26} vs strength  {r:+.3f}   {verdict}")
    print()

    # 3. Convention sensitivity, measured rather than predicted
    print("3. DOES THE MAP SURVIVE A CHANGE OF WEIGHTING CONVENTION?")
    per_conv = {}
    for w in WEIGHTINGS:
        dmg = np.zeros((n_subj, min(args.max_parcels or n_parcel, n_parcel)))
        stg = np.zeros_like(dmg)
        for i, M in enumerate(mats):
            dmg[i], stg[i] = damage_sweep(reweight(M, w), "ac", args.max_parcels)
        per_conv[w] = residualize(dmg, stg, "quadratic").mean(0)
    pairs = [spearmanr(per_conv[a], per_conv[b])[0]
             for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
    stability = float(np.mean(pairs))
    print(f"   degree-corrected map stability across four conventions  {stability:+.3f}")
    if stability > 0.7:
        print("   The ranking is reasonably convention-independent.")
    elif stability > 0.4:
        print("   Partly convention-dependent. Report results under more than one.")
    else:
        print("   CONVENTION-DEPENDENT. A map from any single weighting is not")
        print("   trustworthy on its own here.")
    print()

    # 4. What to do
    print("4. RECOMMENDATION")
    if np.median(tails) >= TAIL_SAFE and stability < 0.7:
        print("   Compress the weight distribution, keeping every edge.")
        log_tail = np.median([tail_ratio(reweight(M, "log")) for M in mats])
        if log_tail < TAIL_SAFE:
            print(f"     log10(1 + w):  tail {np.median(tails):.0f} -> {log_tail:.0f}")
        else:
            print(f"     log10(1 + w) will NOT help here: tail only "
                  f"{np.median(tails):.0f} -> {log_tail:.0f}, because these")
            print(f"     weights are below 1, where log10(1+w) is nearly linear.")
        alpha = recommend_alpha(mats, TAIL_SAFE)
        if alpha is not None:
            achieved = np.median([tail_ratio(power_transform(M, alpha))
                                  for M in mats])
            print(f"     w ** {alpha:.2f}:      tail {np.median(tails):.0f} -> "
                  f"{achieved:.0f}   <- use this")
        else:
            print("     No power transform in (0, 1] reaches the target tail;")
            print("     consider ranking the weights instead.")
        print("   Thresholding also compresses the tail, but discards most of")
        print("   the graph to do it and was measured to work less well.")
    elif stability >= 0.7:
        print("   No action needed on weighting grounds.")
    else:
        print("   Map is unstable but the tail is already compressed, so the")
        print("   cause is elsewhere. Investigate before trusting the ranking.")
    print("\n   Regardless: if section 2 says the measure is degree, the simpler")
    print("   measure is the honest one to report.\n")


if __name__ == "__main__":
    main()
