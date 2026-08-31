"""Does controllability stop being node strength when resections are multi-parcel?

Single-parcel deletion showed d_ave_control correlating with node strength at
+0.909, which makes it a degree proxy and the whole idea redundant. The one
principled reason to expect that to break: strength is additive over a resected
set, while controllability is not. Removing two parcels that are each modestly
connected but jointly carry a bottleneck should cost more controllability than
the sum of their strengths predicts.

This tests that directly, using real HCP-MMP1 surface adjacency so the
resections are anatomically contiguous.
"""

import numpy as np
from scipy.stats import spearmanr

from annot import build_hcp_mmp1_adjacency
from controllability import (
    average_controllability,
    modal_controllability,
    spectral_scale,
)
from data import load_rosen_halgren
from lesion import global_efficiency, grow_resection
from run_real import CSV, partial_spearman


def sweep(A, adjacency, size):
    n = A.shape[0]
    ref = spectral_scale(A)
    intact_ac = average_controllability(A, ref_scale=ref)
    intact_mc = modal_controllability(A, ref_scale=ref)
    baseline_ge = global_efficiency(A)
    strength = A.sum(axis=1)

    rows = {"d_ac": [], "d_mc": [], "d_ge": [], "set_strength": [], "seed": []}
    for seed in range(n):
        resected = grow_resection(A, seed, size, adjacency=adjacency)
        if len(resected) < size:
            continue
        retained = np.setdiff1d(np.arange(n), resected)
        sub = A[np.ix_(retained, retained)]

        rows["seed"].append(seed)
        rows["d_ac"].append(intact_ac[retained].mean()
                            - average_controllability(sub, ref_scale=ref).mean())
        rows["d_mc"].append(intact_mc[retained].mean()
                            - modal_controllability(sub, ref_scale=ref).mean())
        rows["d_ge"].append(baseline_ge - global_efficiency(A, resected))
        # The additive baseline: total connectivity removed.
        rows["set_strength"].append(strength[resected].sum())

    return {k: np.asarray(v) for k, v in rows.items()}


def main():
    A = load_rosen_halgren(CSV)
    adjacency = build_hcp_mmp1_adjacency()

    print("Does controllability escape the degree baseline as resections grow?")
    print("(single-parcel reference: d_ac vs strength = +0.909)\n")
    print(f"{'size':>4}  {'n':>4}  {'dAC~strength':>13}  {'dMC~strength':>13}  "
          f"{'dGE~strength':>13}  {'partial dAC~dGE':>16}")

    for size in (1, 2, 4, 8, 10):
        r = sweep(A, adjacency, size)
        s = r["set_strength"]
        print(f"{size:>4}  {len(s):>4}  "
              f"{spearmanr(r['d_ac'], s)[0]:>+13.3f}  "
              f"{spearmanr(r['d_mc'], s)[0]:>+13.3f}  "
              f"{spearmanr(r['d_ge'], s)[0]:>+13.3f}  "
              f"{partial_spearman(r['d_ac'], r['d_ge'], s):>+16.3f}")


if __name__ == "__main__":
    main()
