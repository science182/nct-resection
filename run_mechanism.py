"""Density or tail heaviness? They are confounded in every result so far.

The density sweep showed that thinning HCP connectomes recovers map stability,
and that was read as density being the cause. But thresholding does two things
at once: it removes edges, lowering density, and it discards the smallest
weights, compressing the weight distribution. Every point in that sweep moved
both knobs together, so the mechanism claim is not identified.

This separates them by building base connectomes on a grid before any weighting
convention is applied:

    alpha   raise every weight to this power. alpha = 1 leaves the original
            heavy tail; smaller alpha compresses it toward uniform without
            removing a single edge, so density is untouched.
    density keep this fraction of edges by magnitude, leaving the surviving
            weights alone.

The four weighting conventions are then applied to each base variant as usual,
and map stability is measured across them. If stability tracks alpha at fixed
density, the tail is the cause and thresholding merely works by compressing it.
If it tracks density at fixed alpha, density is the cause. If both, say so.

Average controllability and strength only; global efficiency is not needed for a
map-stability comparison and costs more than the rest combined.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import itertools  # noqa: E402
import multiprocessing as mp  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from controllability import average_controllability, spectral_scale  # noqa: E402
from data import threshold_edges  # noqa: E402
from run_confounds import residualize  # noqa: E402
from run_weightings import WEIGHTINGS, reweight  # noqa: E402

MAT = "data/streamlineCount.mat"
KEY = "rawStreamlineCounts"
ALPHAS = (1.0, 0.5, 0.25)
DENSITIES = (1.0, 0.30, 0.11)


def base_variant(M, alpha, density):
    """Apply the tail exponent first, then thin, so the two are independent."""
    out = M.copy()
    if alpha != 1.0:
        nz = out > 0
        out[nz] = out[nz] ** alpha
    if density < 1.0:
        out = threshold_edges(out, density)
    return out


def tail_ratio(M):
    """Max over median of the nonzero off-diagonal weights, a crude tail index."""
    off = M[~np.eye(M.shape[0], dtype=bool)]
    nz = off[off > 0]
    return float(nz.max() / np.median(nz)) if nz.size else float("nan")


def ac_sweep(A):
    """Per-parcel deletion damage in average controllability, plus strength."""
    n = A.shape[0]
    ref = spectral_scale(A)
    intact = average_controllability(A, ref_scale=ref)
    d_ac = np.zeros(n)
    for i in range(n):
        keep = np.setdiff1d(np.arange(n), [i])
        d_ac[i] = intact[keep].mean() - average_controllability(
            A[np.ix_(keep, keep)], ref_scale=ref).mean()
    return d_ac, A.sum(axis=1)


def _one(args):
    index, alpha, density = args
    with h5py.File(MAT, "r") as f:
        M = np.asarray(f[KEY][index], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    M = (M + M.T) / 2.0
    base = base_variant(M, alpha, density)
    out = {}
    for w in WEIGHTINGS:
        out[w] = ac_sweep(reweight(base, w))
    return index, alpha, density, out, tail_ratio(base), float(np.mean(
        base[~np.eye(360, dtype=bool)] > 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()

    cells = list(itertools.product(ALPHAS, DENSITIES))
    tasks = [(i, a, d) for a, d in cells for i in range(args.n)]
    store = {(a, d): {f"{w}_{k}": np.zeros((args.n, 360))
                      for w in WEIGHTINGS for k in ("ac", "st")}
             for a, d in cells}
    meta = {}

    print(f"{len(cells)} cells x {args.n} subjects, {args.workers} workers\n",
          flush=True)
    with mp.Pool(args.workers) as pool:
        for done, (index, alpha, density, res, tr, dens) in enumerate(
                pool.imap_unordered(_one, tasks), 1):
            for w, (ac, st) in res.items():
                store[(alpha, density)][f"{w}_ac"][index] = ac
                store[(alpha, density)][f"{w}_st"][index] = st
            meta[(alpha, density)] = (tr, dens)
            if done % 12 == 0 or done == len(tasks):
                print(f"  {done}/{len(tasks)}", flush=True)

    print("\nmap stability across the four weighting conventions")
    print(f"  {'alpha':>6} {'density set':>12} {'actual':>8} {'tail max/med':>13} "
          f"{'stability':>11}")
    table = {}
    for alpha, density in cells:
        s = store[(alpha, density)]
        maps = {w: residualize(s[f"{w}_ac"], s[f"{w}_st"], "quadratic").mean(0)
                for w in WEIGHTINGS}
        vals = [spearmanr(maps[a], maps[b])[0]
                for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
        table[(alpha, density)] = float(np.mean(vals))
        tr, dens = meta[(alpha, density)]
        print(f"  {alpha:>6.2f} {density:>12.2f} {dens:>8.2f} {tr:>13.0f} "
              f"{np.mean(vals):>+11.3f}")

    print("\n  effect of compressing the tail, at fixed density:")
    for d in DENSITIES:
        deltas = table[(ALPHAS[-1], d)] - table[(ALPHAS[0], d)]
        print(f"    density {d:.2f}: alpha 1.00 -> 0.25 changes stability by "
              f"{deltas:+.3f}")
    print("\n  effect of thinning, at fixed tail:")
    for a in ALPHAS:
        deltas = table[(a, DENSITIES[-1])] - table[(a, DENSITIES[0])]
        print(f"    alpha {a:.2f}: density 1.00 -> 0.11 changes stability by "
              f"{deltas:+.3f}")


if __name__ == "__main__":
    main()
