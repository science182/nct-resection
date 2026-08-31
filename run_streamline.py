"""Replicate the language finding on raw streamline counts.

Everything so far used Fpt, the fraction of streamlines reaching a target. Lin
et al. 2024 used raw streamline counts instead, so this checks whether the
result depends on that choice.

The two are not simply rescaled versions of each other. Fpt is normalized per
seed parcel, which removes the effect of a parcel's size on how many streamlines
it emits. Raw counts keep it: strength spans roughly twentyfold across parcels
here against about 2.4-fold for Fpt. So raw counts carry a stronger built-in
size confound, and this is a genuinely different test rather than a reweighting.

Raw counts are also asymmetric, since seeding from i and finding j is not the
same measurement as the reverse, so they are symmetrized before use.

Running 48 subjects rather than 24 allows an independent split-half
replication, which is the weakest point in the Fpt result.
"""

import os
import sys

import h5py
import numpy as np

from run_subject_deletions import sweep_subject

MAT = "data/streamlineCount.mat"
KEY = "rawStreamlineCounts"
OUT = "data/streamline_deletions.npz"


def load_subject(dset, i):
    M = np.asarray(dset[i], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    return (M + M.T) / 2.0  # raw counts are directional


def main():
    """Checkpoints after every subject.

    An earlier version only wrote at the end, and the process was killed part
    way through, losing the whole run. `done` records how many rows are real so
    an interrupted run can be resumed or analysed as far as it got.
    """
    n_subj = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    with h5py.File(MAT, "r") as f:
        dset = f[KEY]
        n_parcel = dset.shape[1]

        D_AC = np.zeros((n_subj, n_parcel))
        D_GE = np.zeros((n_subj, n_parcel))
        ST = np.zeros((n_subj, n_parcel))
        start = 0

        if os.path.exists(OUT):
            prev = np.load(OUT)
            done = int(prev["done"])
            if done and prev["d_ac"].shape[1] == n_parcel:
                k = min(done, n_subj)
                D_AC[:k], D_GE[:k], ST[:k] = (prev["d_ac"][:k], prev["d_ge"][:k],
                                              prev["strength"][:k])
                start = k
                print(f"resuming from subject {start}", flush=True)

        for s in range(start, n_subj):
            D_AC[s], D_GE[s], ST[s] = sweep_subject(load_subject(dset, s))
            np.savez_compressed(OUT, d_ac=D_AC, d_ge=D_GE, strength=ST,
                                done=s + 1)
            print(f"  {s + 1}/{n_subj}", flush=True)

    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
