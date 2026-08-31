"""Does controllability carry individual information that node strength does not?

Everything before this ran on a group average, which by construction cannot
speak to interindividual variability. That matters because the finding in Lin
et al. 2024 is precisely about individual differences: the parcel whose removal
costs the most global efficiency is not the same parcel in every subject.

Average controllability turned out to track node strength at +0.909 on the
group average. That does not settle the individual question. A measure can be a
near-deterministic function of degree on the mean brain and still deviate from
degree in a subject-specific, reproducible way.

Three tests, in increasing order of how much they would prove:

1. Within-subject redundancy. Spearman(nodal AC, nodal strength) per subject.
   If this sits near 1.0 for everyone, AC is degree for everyone.

2. Deviation structure. Regress each parcel's AC on its strength across
   subjects. The residual is AC variation strength cannot explain. If residuals
   are pure noise, there is nothing individual here.

3. Fingerprinting, the real test. Split parcels into two halves. Correlate each
   subject's residual profile on half A against every subject's on half B. If
   the self-match wins, the residual is a stable individual signature rather
   than noise, which is what "connectotype" has to mean to be useful.

Full deletion sweeps are not feasible here: one subject takes about two
minutes, so 1065 would take a day and a half. Nodal controllability is one
decomposition per subject.

Read the fingerprinting result carefully. This dataset has one connectome per
subject, so the split is across parcels within a single scan, not across
sessions. Beating chance therefore shows the residual is spatially structured
and subject-specific in that scan. It does NOT show the residual is a stable
trait, because subject-specific tractography artifacts would produce the same
result. Test-retest data is what separates those, and this cannot substitute.
A failure here is still decisive in the other direction: a residual that cannot
be fingerprinted even within a single scan is noise.
"""

import sys

import h5py
import numpy as np
from scipy.stats import spearmanr

from controllability import average_controllability, modal_controllability

MAT = "data/individualConnectivity.mat"
CACHE = "data/subject_profiles.npz"


def inspect(path=MAT):
    with h5py.File(path, "r") as f:
        keys = [k for k in f.keys() if not k.startswith("#")]
        for k in keys:
            print(f"  {k}: shape={f[k].shape} dtype={f[k].dtype}")
    return keys


def subject_slice(dset, i):
    """Read one subject's 360x360 matrix, whichever axis subjects live on."""
    shape = dset.shape
    axis = int(np.argmax([s for s in shape]))  # subjects axis is the long one
    if axis == 0:
        M = np.asarray(dset[i], dtype=float)
    elif axis == 1:
        M = np.asarray(dset[:, i], dtype=float)
    else:
        M = np.asarray(dset[:, :, i], dtype=float)
    return M


def build_profiles(path=MAT, cache=CACHE, limit=None):
    """Compute nodal AC, MC and strength for every subject."""
    with h5py.File(path, "r") as f:
        key = [k for k in f.keys() if not k.startswith("#")][0]
        dset = f[key]
        shape = dset.shape
        n_subj = max(shape)
        n_parcel = min(shape)
        if limit:
            n_subj = min(n_subj, limit)
        print(f"  {key}: {shape} -> {n_subj} subjects, {n_parcel} parcels")

        ac = np.zeros((n_subj, n_parcel))
        mc = np.zeros((n_subj, n_parcel))
        st = np.zeros((n_subj, n_parcel))

        for i in range(n_subj):
            M = subject_slice(dset, i)
            if M.shape != (n_parcel, n_parcel):
                M = M.reshape(n_parcel, n_parcel)
            M[~np.isfinite(M)] = 0.0
            np.fill_diagonal(M, 0.0)
            M = (M + M.T) / 2.0
            if not np.any(M):
                continue
            ac[i] = average_controllability(M)
            mc[i] = modal_controllability(M)
            st[i] = M.sum(axis=1)
            if (i + 1) % 100 == 0:
                print(f"    {i + 1}/{n_subj}", flush=True)

    np.savez_compressed(cache, ac=ac, mc=mc, strength=st)
    return ac, mc, st


def residualize(profile, strength):
    """Remove, parcel by parcel across subjects, whatever strength explains."""
    resid = np.zeros_like(profile)
    for p in range(profile.shape[1]):
        x = strength[:, p]
        y = profile[:, p]
        design = np.column_stack([np.ones_like(x), x])
        beta = np.linalg.lstsq(design, y, rcond=None)[0]
        resid[:, p] = y - design @ beta
    return resid


def fingerprint(resid, seed=0):
    """Can a subject be identified from their residual profile?

    Split parcels in half, match each subject's half-A profile against every
    subject's half-B profile, and count how often the correct subject wins.
    Chance is 1/n_subjects.
    """
    rng = np.random.default_rng(seed)
    n_subj, n_parcel = resid.shape
    perm = rng.permutation(n_parcel)
    a, b = perm[: n_parcel // 2], perm[n_parcel // 2:]

    A = resid[:, a]
    B = resid[:, b]
    A = (A - A.mean(1, keepdims=True)) / (A.std(1, keepdims=True) + 1e-12)
    B = (B - B.mean(1, keepdims=True)) / (B.std(1, keepdims=True) + 1e-12)
    sim = A @ B.T / A.shape[1]

    best = np.argmax(sim, axis=1)
    accuracy = float((best == np.arange(n_subj)).mean())
    self_sim = np.diag(sim)
    off = sim[~np.eye(n_subj, dtype=bool)]
    return accuracy, float(self_sim.mean()), float(off.mean())


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    print("file contents:")
    inspect()
    print("\nbuilding per-subject profiles...")
    ac, mc, st = build_profiles(limit=limit)
    n = ac.shape[0]

    print(f"\n=== 1. within-subject redundancy (n={n}) ===")
    rho = np.array([spearmanr(ac[i], st[i])[0] for i in range(n)])
    rho_mc = np.array([spearmanr(mc[i], st[i])[0] for i in range(n)])
    print(f"  Spearman(AC, strength): mean {rho.mean():+.3f}, "
          f"sd {rho.std():.3f}, range [{rho.min():+.3f}, {rho.max():+.3f}]")
    print(f"  Spearman(MC, strength): mean {rho_mc.mean():+.3f}, "
          f"sd {rho_mc.std():.3f}")

    print("\n=== 2. deviation structure ===")
    resid_ac = residualize(ac, st)
    frac = 1.0 - resid_ac.var(axis=0).sum() / ac.var(axis=0).sum()
    print(f"  variance in AC explained by strength: {frac:.1%}")
    print(f"  residual sd (mean over parcels): {resid_ac.std(axis=0).mean():.4g}")

    print("\n=== 3. fingerprinting the residual ===")
    for name, prof in (("AC residual (strength removed)", resid_ac),
                       ("AC raw", ac),
                       ("strength", st)):
        acc, s_self, s_off = fingerprint(prof)
        print(f"  {name:<32} accuracy {acc:6.1%}  "
              f"(chance {1 / n:.2%})  self r={s_self:+.3f} other r={s_off:+.3f}")


if __name__ == "__main__":
    main()
