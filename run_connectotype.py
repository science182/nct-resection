"""Does controllability capture individual variation that node strength does not?

Replaces an invalid test. The first attempt split parcels into halves and
correlated a subject's profile on half A against profiles on half B. Those
vectors are indexed by different parcels, so they have no reason to align even
within a subject, and the test returned chance accuracy for every input
including node strength. Strength profiles are reliably identifiable in the
fingerprinting literature, so a design that cannot recover them is broken.

True fingerprinting needs two sessions per subject. This dataset has one
connectome per subject, so that is simply not available. These tests are what
the data can actually support:

A. Similarity structure. Build the subject-by-subject correlation matrix from
   controllability profiles and from strength profiles. If the two carry the
   same information about who resembles whom, they are interchangeable.

B. Connectotype diversity. Lin et al. 2024 found that the most damaging parcel
   differs between individuals. Count how many distinct parcels take the top
   rank across subjects, by controllability and by strength.

C. Structure versus noise. Compare the residual's between-subject similarity
   against a null built by shuffling parcels within each subject. Real
   individual structure survives; noise does not.

D. The nonlinearity confound. Residualizing on strength linearly leaves any
   nonlinear dependence behind, so a residual could re-encode degree and look
   novel. Check for leftover monotonic dependence directly.
"""

import numpy as np
from scipy.stats import spearmanr

CACHE = "data/subject_profiles.npz"


def zscore_rows(X):
    return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)


def similarity(profiles):
    """Subject-by-subject correlation across the parcel profile."""
    Z = zscore_rows(profiles)
    return Z @ Z.T / Z.shape[1]


def offdiag(M):
    return M[~np.eye(M.shape[0], dtype=bool)]


def residualize(profile, strength):
    resid = np.zeros_like(profile)
    for p in range(profile.shape[1]):
        x, y = strength[:, p], profile[:, p]
        design = np.column_stack([np.ones_like(x), x])
        resid[:, p] = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return resid


def main():
    d = np.load(CACHE)
    ac, mc, st = d["ac"], d["mc"], d["strength"]
    n_subj, n_parcel = ac.shape
    print(f"{n_subj} subjects, {n_parcel} parcels\n")

    resid = residualize(ac, st)

    print("=== D. does the residual still encode strength nonlinearly? ===")
    rhos = np.array([spearmanr(resid[:, p], st[:, p])[0] for p in range(n_parcel)])
    print(f"  Spearman(residual, strength) per parcel: mean {rhos.mean():+.3f}, "
          f"max |rho| {np.abs(rhos).max():.3f}")
    print(f"  parcels with |rho| > 0.2: {(np.abs(rhos) > 0.2).sum()} / {n_parcel}")
    verdict = "clean" if np.abs(rhos).max() < 0.2 else "LEAKY, residual re-encodes degree"
    print(f"  -> {verdict}\n")

    print("=== A. do controllability and strength say the same thing about who is similar? ===")
    sim_ac = similarity(ac)
    sim_st = similarity(st)
    sim_rs = similarity(resid)
    r_raw = spearmanr(offdiag(sim_ac), offdiag(sim_st))[0]
    r_res = spearmanr(offdiag(sim_rs), offdiag(sim_st))[0]
    print(f"  subject-similarity agreement, AC vs strength:        {r_raw:+.3f}")
    print(f"  subject-similarity agreement, AC residual vs strength:{r_res:+.3f}")
    print(f"  -> {'redundant' if abs(r_raw) > 0.9 else 'AC similarity is not strength similarity'}\n")

    print("=== B. connectotype diversity: which parcel ranks top? ===")
    for name, prof in (("average controllability", ac), ("strength", st),
                       ("AC residual", resid)):
        top = np.argmax(prof, axis=1)
        counts = np.bincount(top, minlength=n_parcel)
        distinct = int((counts > 0).sum())
        modal = counts.max() / n_subj
        print(f"  {name:<24} {distinct:>4} distinct top parcels, "
              f"most common holds {modal:.1%} of subjects")
    print()

    print("=== C. is the residual structured, or noise? ===")
    rng = np.random.default_rng(0)
    shuffled = np.array([rng.permutation(row) for row in resid])
    real_off = offdiag(similarity(resid))
    null_off = offdiag(similarity(shuffled))
    print(f"  real residual similarity:    sd {real_off.std():.4f}, "
          f"|r|>0.3 in {(np.abs(real_off) > 0.3).mean():.2%} of pairs")
    print(f"  parcel-shuffled null:        sd {null_off.std():.4f}, "
          f"|r|>0.3 in {(np.abs(null_off) > 0.3).mean():.2%} of pairs")
    ratio = real_off.std() / null_off.std()
    print(f"  spread ratio real/null: {ratio:.2f}")
    print(f"  -> {'structured beyond chance' if ratio > 1.5 else 'indistinguishable from noise'}")


if __name__ == "__main__":
    main()
