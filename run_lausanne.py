"""Independent replication in a second dataset.

Everything else here uses Rosen & Halgren HCP connectomes, which makes the
weighting result a statement about one tractography pipeline rather than about
the method. This repeats the core test in a dataset that differs on every axis
that matters:

    subjects      70 healthy adults, not the HCP 1065
    site          different acquisition
    pipeline      Lausanne / Connectome Mapper, not probtrackx
    parcellation  Lausanne 219, not HCP-MMP1 360
    edge weight   fiber density, normalized by streamline length and region
                  surface area, not fractional probability or raw counts
    density       0.111, sparse, against 1.000 fully dense for HCP

Zenodo 2872624, CC-BY 4.0.

The dataset ships one native weighting, so the four conventions are derived from
it exactly as run_weightings.py does for streamline counts. The question is
whether the degree-corrected damage map is as convention-dependent here as it is
in HCP. Network labels are not available for the Lausanne atlas, so this tests
map stability rather than network enrichment, which is the actual claim anyway.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
import scipy.io as sio  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from run_confounds import residualize  # noqa: E402
from run_subject_deletions import sweep_subject  # noqa: E402
from run_weightings import WEIGHTINGS, reweight  # noqa: E402

MAT = "data/lausanne70.mat"
SCALES = {68: 0, 114: 1, 219: 2, 448: 3, 1000: 4}
OUT = "data/lausanne_deletions.npz"


def load_scale(scale):
    m = sio.loadmat(MAT)
    A = np.asarray(m["connMatrices"][0, 0]["SC"][SCALES[scale], 0], dtype=float)
    return A  # (parcels, parcels, subjects)


def _one(args):
    flat, n_parcel, index = args
    M = flat.reshape(n_parcel, n_parcel)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    M = (M + M.T) / 2.0
    out = {}
    for w in WEIGHTINGS:
        ac, ge, st = sweep_subject(reweight(M, w))
        out[w] = (ac, ge, st)
    return index, out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", type=int, default=219, choices=sorted(SCALES))
    ap.add_argument("--n", type=int, default=70)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()

    A = load_scale(args.scale)
    n_parcel, n_subj = A.shape[0], min(args.n, A.shape[2])
    print(f"Lausanne scale {args.scale}: {n_parcel} parcels, {n_subj} subjects")

    store = {f"{w}_{k}": np.zeros((n_subj, n_parcel))
             for w in WEIGHTINGS for k in ("d_ac", "d_ge", "strength")}

    tasks = [(A[:, :, i].ravel().copy(), n_parcel, i) for i in range(n_subj)]
    with mp.Pool(args.workers) as pool:
        for done, (index, res) in enumerate(pool.imap_unordered(_one, tasks), 1):
            for w, (ac, ge, st) in res.items():
                store[f"{w}_d_ac"][index] = ac
                store[f"{w}_d_ge"][index] = ge
                store[f"{w}_strength"][index] = st
            if done % 10 == 0 or done == n_subj:
                print(f"  {done}/{n_subj}", flush=True)

    np.savez_compressed(OUT, scale=args.scale, n_subj=n_subj, **store)

    print(f"\n=== degree confound, Lausanne {args.scale} ===")
    for w in WEIGHTINGS:
        ac, st = store[f"{w}_d_ac"], store[f"{w}_strength"]
        r = np.mean([spearmanr(ac[s], st[s])[0] for s in range(n_subj)])
        print(f"  {w:<9} d_ac vs strength within subject: {r:+.3f}")

    print(f"\n=== map stability across weighting conventions ===")
    maps = {w: residualize(store[f"{w}_d_ac"], store[f"{w}_strength"],
                           "quadratic") for w in WEIGHTINGS}
    raw = {w: store[f"{w}_d_ac"] for w in WEIGHTINGS}
    for label, group in (("raw damage", raw), ("degree-corrected", maps)):
        vals = [spearmanr(group[a].mean(0), group[b].mean(0))[0]
                for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
        print(f"  {label:<18} mean pairwise spatial agreement: {np.mean(vals):+.3f}")
    strength = {w: store[f"{w}_strength"] for w in WEIGHTINGS}
    vals = [spearmanr(strength[a].mean(0), strength[b].mean(0))[0]
            for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
    print(f"  {'node strength':<18} mean pairwise spatial agreement: {np.mean(vals):+.3f}")
    print(f"\n  (HCP for comparison: degree-corrected +0.289, strength +0.748)")


if __name__ == "__main__":
    main()
