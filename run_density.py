"""Does connectome density explain why the weighting instability did not replicate?

In HCP the degree-corrected map collapsed across weighting conventions (+0.289).
In the Lausanne dataset it held up (+0.756). The two differ in many ways, but the
most obvious is structural: HCP probtrackx output is fully dense with edge
weights spanning about 250,000 to 1, while Lausanne is 11 percent dense spanning
about 20,000 to 1.

If density and tail heaviness are the cause, then thinning HCP down toward
Lausanne's density should make its map stability climb toward Lausanne's number.
If stability stays low, density is not the explanation and the difference lies
somewhere else in the pipeline.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from data import threshold_edges  # noqa: E402
from run_confounds import residualize  # noqa: E402
from run_subject_deletions import sweep_subject  # noqa: E402
from run_weightings import WEIGHTINGS, reweight  # noqa: E402

MAT = "data/streamlineCount.mat"
KEY = "rawStreamlineCounts"


def _one(args):
    index, keep = args
    with h5py.File(MAT, "r") as f:
        M = np.asarray(f[KEY][index], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    M = (M + M.T) / 2.0
    if keep < 1.0:
        M = threshold_edges(M, keep)
    out = {}
    for w in WEIGHTINGS:
        ac, _, st = sweep_subject(reweight(M, w))
        out[w] = (ac, st)
    return index, out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()

    print("HCP streamline counts, thinned toward the Lausanne density\n")
    print(f"  {'kept':>6} {'density':>9} {'corrected map stability':>25}")

    for keep in (1.0, 0.50, 0.30, 0.11):
        store = {f"{w}_{k}": np.zeros((args.n, 360))
                 for w in WEIGHTINGS for k in ("d_ac", "strength")}
        tasks = [(i, keep) for i in range(args.n)]
        with mp.Pool(args.workers) as pool:
            for index, res in pool.imap_unordered(_one, tasks):
                for w, (ac, st) in res.items():
                    store[f"{w}_d_ac"][index] = ac
                    store[f"{w}_strength"][index] = st

        maps = {w: residualize(store[f"{w}_d_ac"], store[f"{w}_strength"],
                               "quadratic").mean(0) for w in WEIGHTINGS}
        vals = [spearmanr(maps[a], maps[b])[0]
                for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]

        with h5py.File(MAT, "r") as f:
            probe = np.asarray(f[KEY][0], dtype=float)
        probe = (probe + probe.T) / 2.0
        np.fill_diagonal(probe, 0.0)
        if keep < 1.0:
            probe = threshold_edges(probe, keep)
        off = probe[~np.eye(360, dtype=bool)]
        print(f"  {keep:>6.0%} {np.mean(off > 0):>9.2f} {np.mean(vals):>+25.3f}",
              flush=True)

    print(f"\n  Lausanne 219 at density 0.11 for comparison: +0.756")
    print(f"  n = {args.n} subjects per condition")


if __name__ == "__main__":
    main()
