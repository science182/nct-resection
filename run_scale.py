"""Per-subject deletion sweeps at full scale.

The original sweep did 24 subjects in 9.4 hours, which is why every per-subject
result in this project rests on 24 of the 1065 available. Three changes fix
that, and none of them alter the numbers (verified identical to 3e-12 against
the original run):

  1. The connectome is symmetric, so `eigh` replaces the general `schur`.
     About six times faster.
  2. The graph is dense, so Floyd-Warshall replaces Dijkstra for all-pairs
     shortest paths. About five times faster.
  3. BLAS threading is switched off and parallelism moved to whole subjects.
     At 359x359 the matrices are small enough that threading costs more than it
     returns: one thread measured 92 ms per deletion against 140 ms on eight.
     One thread per worker, one subject per worker, scales close to linearly.

Together that is roughly 300x, which turns "24 subjects overnight" into "all
1065 in about an hour".

    python3 run_scale.py --n 1065 --workers 8
    python3 run_scale.py --n 200 --data streamline
"""

import os

# Must precede the numpy import so BLAS picks it up in every worker.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import time  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402

from run_subject_deletions import sweep_subject  # noqa: E402

SOURCES = {
    "fpt": ("data/individualConnectivity.mat", "individualConnectivity"),
    "streamline": ("data/streamlineCount.mat", "rawStreamlineCounts"),
}


def _load(path, key, index):
    """Open per call. An h5py handle cannot be shared across processes."""
    with h5py.File(path, "r") as f:
        M = np.asarray(f[key][index], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    return (M + M.T) / 2.0  # raw streamline counts are directional


def _one(args):
    path, key, index = args
    ac, ge, st = sweep_subject(_load(path, key, index))
    return index, ac, ge, st


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=1065)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--data", choices=sorted(SOURCES), default="fpt")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    path, key = SOURCES[args.data]
    out = args.out or f"data/scale_{args.data}_deletions.npz"

    with h5py.File(path, "r") as f:
        available, n_parcel = f[key].shape[0], f[key].shape[1]
    n = min(args.n, available)

    AC = np.zeros((n, n_parcel))
    GE = np.zeros((n, n_parcel))
    ST = np.zeros((n, n_parcel))
    done = np.zeros(n, dtype=bool)

    if os.path.exists(out):
        prev = np.load(out)
        if prev["d_ac"].shape == AC.shape:
            AC, GE, ST = prev["d_ac"], prev["d_ge"], prev["strength"]
            done = prev["done_mask"]
            print(f"resuming, {int(done.sum())}/{n} already complete", flush=True)

    todo = [(path, key, i) for i in range(n) if not done[i]]
    if not todo:
        print("nothing to do")
        return

    print(f"{len(todo)} subjects on {args.workers} workers "
          f"({args.data}, {n_parcel} parcels)", flush=True)
    start = time.time()
    completed = 0

    with mp.Pool(args.workers) as pool:
        for index, ac, ge, st in pool.imap_unordered(_one, todo, chunksize=1):
            AC[index], GE[index], ST[index] = ac, ge, st
            done[index] = True
            completed += 1
            if completed % 10 == 0 or completed == len(todo):
                np.savez_compressed(out, d_ac=AC, d_ge=GE, strength=ST,
                                    done_mask=done, done=int(done.sum()))
                rate = completed / (time.time() - start)
                left = (len(todo) - completed) / rate if rate else 0
                print(f"  {completed}/{len(todo)}  "
                      f"{rate * 60:.1f}/min  ~{left / 60:.0f} min left", flush=True)

    np.savez_compressed(out, d_ac=AC, d_ge=GE, strength=ST,
                        done_mask=done, done=int(done.sum()))
    print(f"saved {out} in {(time.time() - start) / 60:.1f} min")


if __name__ == "__main__":
    main()
