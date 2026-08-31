"""Per-subject deletion sweeps: does resection risk itself vary individually?

The connectotype analysis showed nodal controllability carries subject-specific
structure beyond node strength. That is not yet the claim that matters. The
claim that matters is about resection damage: whether the parcel whose removal
costs the most controllability differs between individuals, and whether that
variation is something global efficiency and degree do not already predict.

That requires deleting every parcel in every subject, which costs about 90
seconds per subject, so this runs on a subset rather than all 1065.
"""

import sys

import h5py
import numpy as np
from scipy.stats import spearmanr

from controllability import average_controllability, spectral_scale
from lesion import global_efficiency

MAT = "data/individualConnectivity.mat"
OUT = "data/subject_deletions.npz"


def load_subject(dset, i):
    M = np.asarray(dset[i], dtype=float)
    M[~np.isfinite(M)] = 0.0
    np.fill_diagonal(M, 0.0)
    return (M + M.T) / 2.0


def sweep_subject(A):
    n = A.shape[0]
    ref = spectral_scale(A)
    intact_ac = average_controllability(A, ref_scale=ref)
    baseline_ge = global_efficiency(A)

    d_ac = np.zeros(n)
    d_ge = np.zeros(n)
    for i in range(n):
        retained = np.setdiff1d(np.arange(n), [i])
        sub = A[np.ix_(retained, retained)]
        d_ac[i] = intact_ac[retained].mean() - average_controllability(
            sub, ref_scale=ref).mean()
        d_ge[i] = baseline_ge - global_efficiency(A, [i])
    return d_ac, d_ge, A.sum(axis=1)


def main():
    n_subj = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    with h5py.File(MAT, "r") as f:
        dset = f["individualConnectivity"]
        n_parcel = dset.shape[1]
        D_AC = np.zeros((n_subj, n_parcel))
        D_GE = np.zeros((n_subj, n_parcel))
        ST = np.zeros((n_subj, n_parcel))
        for s in range(n_subj):
            D_AC[s], D_GE[s], ST[s] = sweep_subject(load_subject(dset, s))
            print(f"  subject {s + 1}/{n_subj}", flush=True)

    np.savez_compressed(OUT, d_ac=D_AC, d_ge=D_GE, strength=ST)

    print(f"\n=== per-subject deletion results, n={n_subj} ===")
    rho = np.array([spearmanr(D_AC[s], ST[s])[0] for s in range(n_subj)])
    rho_ge = np.array([spearmanr(D_AC[s], D_GE[s])[0] for s in range(n_subj)])
    print(f"  Spearman(d_ac, strength) within subject: mean {rho.mean():+.3f} "
          f"sd {rho.std():.3f}")
    print(f"  Spearman(d_ac, d_ge)     within subject: mean {rho_ge.mean():+.3f} "
          f"sd {rho_ge.std():.3f}")

    print("\n  which parcel is worst, across subjects:")
    for name, M in (("d_ave_control", D_AC), ("d_global_eff", D_GE),
                    ("strength", ST)):
        top = np.argmax(M, axis=1)
        counts = np.bincount(top, minlength=n_parcel)
        print(f"    {name:<15} {int((counts > 0).sum()):>3} distinct winners, "
              f"modal parcel holds {counts.max() / n_subj:.0%}")

    # Is a subject's controllability ranking better predicted by their own GE
    # ranking, or by another subject's? If own > others, the individual
    # deviation is real rather than shared structure.
    own, other = [], []
    for s in range(n_subj):
        own.append(spearmanr(D_AC[s], D_GE[s])[0])
        for t in range(n_subj):
            if t != s:
                other.append(spearmanr(D_AC[s], D_GE[t])[0])
    print(f"\n  d_ac vs own-subject d_ge:   {np.mean(own):+.4f}")
    print(f"  d_ac vs other-subject d_ge: {np.mean(other):+.4f}")
    print(f"  individual-specific gain:   {np.mean(own) - np.mean(other):+.4f}")


if __name__ == "__main__":
    main()
