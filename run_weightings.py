"""Which measures survive a change of edge-weighting convention?

Fpt weighting put the degree-corrected controllability peak on language cortex.
Raw streamline counts put it on ventromedial visual cortex, with the two maps
uncorrelated. Since Fpt is defined as streamline counts normalized per seed
parcel, the difference between the two has to come from that normalization.

This derives four weightings from the same streamline count matrices, so the
underlying tractography is held fixed and only the convention changes:

    raw          counts as stored
    log          log10(1 + counts), compresses the ~20x strength range
    rownorm      counts divided by row sum, which is Fpt reconstructed
    binary       top 30 percent of edges, weight discarded entirely

Then it asks two questions. Does the degree-corrected controllability map agree
across conventions? And does global efficiency, the measure this project set out
to improve on, hold up any better under the same test?

A measure whose map depends on the convention cannot support a clinical score,
whatever its p-values look like.
"""

import os
import sys

import h5py
import numpy as np

from run_subject_deletions import sweep_subject

MAT = "data/streamlineCount.mat"
KEY = "rawStreamlineCounts"
OUT = "data/weighting_deletions.npz"
WEIGHTINGS = ("raw", "log", "rownorm", "binary")


def reweight(M, kind):
    """M arrives symmetrized with a zero diagonal."""
    if kind == "raw":
        return M
    if kind == "log":
        return np.log10(1.0 + M)
    if kind == "rownorm":
        rs = M.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        R = M / rs
        return (R + R.T) / 2.0
    if kind == "binary":
        iu = np.triu_indices_from(M, k=1)
        w = M[iu]
        nz = w[w > 0]
        cut = np.quantile(nz, 0.70)
        B = (M >= cut).astype(float)
        np.fill_diagonal(B, 0.0)
        return B
    raise ValueError(kind)


def load_subject(dset, i):
    M = np.asarray(dset[i], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    return (M + M.T) / 2.0


def main():
    n_subj = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    store = {}
    start = 0
    if os.path.exists(OUT):
        prev = np.load(OUT)
        if int(prev["done"]) and int(prev["n_subj"]) == n_subj:
            store = {k: prev[k] for k in prev.files
                     if k not in ("done", "n_subj")}
            start = int(prev["done"])
            print(f"resuming from subject {start}", flush=True)

    for w in WEIGHTINGS:
        for field in ("d_ac", "d_ge", "strength"):
            store.setdefault(f"{w}_{field}", np.zeros((n_subj, 360)))

    with h5py.File(MAT, "r") as f:
        dset = f[KEY]
        for s in range(start, n_subj):
            base = load_subject(dset, s)
            for w in WEIGHTINGS:
                ac, ge, st = sweep_subject(reweight(base, w))
                store[f"{w}_d_ac"][s] = ac
                store[f"{w}_d_ge"][s] = ge
                store[f"{w}_strength"][s] = st
            np.savez_compressed(OUT, done=s + 1, n_subj=n_subj, **store)
            print(f"  {s + 1}/{n_subj}", flush=True)

    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
