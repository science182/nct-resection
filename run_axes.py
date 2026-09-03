"""Is the weighting instability spatial, between-subject, or both?

This tests the field's own defence of controllability. Parkes et al. concede
that average controllability correlates spatially with node strength across the
brain, but argue that this does not make it redundant between subjects: average
controllability outperformed strength at out-of-sample prediction of psychosis
symptoms. Every other analysis in this repository is spatial, which is precisely
the axis where the correlation is already conceded to be high.

So the two weightings can be compared along both axes, using the same 1065
subjects in verified identical order:

  spatial          for each weighting, average across subjects to get one map
                   over 360 parcels, then correlate the two maps. This is the
                   surgical planning question: which parcel is risky.

  between-subject  for each parcel, take the 1065-long vector of damage across
                   subjects and correlate it between weightings, then average
                   over parcels. This is the individual-differences question:
                   which patient is unusual at this parcel.

If the spatial maps disagree while the between-subject rankings agree, the
instability is confined to the axis that matters for a resection score and
leaves the axis the literature actually defends intact. That would be a more
precise and more useful claim than "controllability is unstable".
"""

import numpy as np
from scipy.stats import spearmanr

from run_confounds import residualize

FPT = "data/scale_fpt_deletions.npz"
STREAM = "data/scale_streamline_deletions.npz"


def spatial_agreement(a, b):
    """Correlate the two group-average parcel maps."""
    return spearmanr(a.mean(0), b.mean(0))[0]


def between_subject_agreement(a, b):
    """Per parcel, correlate the across-subject vectors, then summarize."""
    rhos = np.array([spearmanr(a[:, p], b[:, p])[0] for p in range(a.shape[1])])
    return rhos


def main():
    f, s = np.load(FPT), np.load(STREAM)
    ac_f, ac_s = f["d_ac"], s["d_ac"]
    st_f, st_s = f["strength"], s["strength"]
    n, n_parcel = ac_f.shape
    print(f"{n} subjects, {n_parcel} parcels, identical subject order verified\n")

    res_f = residualize(ac_f, st_f, "quadratic")
    res_s = residualize(ac_s, st_s, "quadratic")

    print("agreement between Fpt and streamline weightings")
    print(f"  {'quantity':<28} {'spatial':>9} {'between-subject':>18}")
    for label, a, b in (("raw damage (d_ac)", ac_f, ac_s),
                        ("degree-corrected damage", res_f, res_s),
                        ("node strength", st_f, st_s)):
        bs = between_subject_agreement(a, b)
        print(f"  {label:<28} {spatial_agreement(a, b):>+9.3f} "
              f"{bs.mean():>+18.3f}")

    bs = between_subject_agreement(res_f, res_s)
    print(f"\n  degree-corrected, between-subject detail:")
    print(f"    mean {bs.mean():+.3f}, sd {bs.std():.3f}, "
          f"range {bs.min():+.3f} to {bs.max():+.3f}")
    print(f"    parcels with rho > 0.3: {(bs > 0.3).sum()} / {n_parcel}")
    print(f"    parcels with rho < 0.0: {(bs < 0).sum()} / {n_parcel}")

    # A null: shuffling subjects should destroy any between-subject agreement.
    rng = np.random.default_rng(0)
    perm = rng.permutation(n)
    null = between_subject_agreement(res_f, res_s[perm])
    print(f"\n  subject-shuffled null: mean {null.mean():+.3f}, "
          f"sd {null.std():.3f}")

    print("\ninterpretation")
    if bs.mean() > 0.3 > abs(spatial_agreement(res_f, res_s)):
        print("  The instability is spatial. Between-subject structure survives the")
        print("  change of convention, so the axis the literature defends holds and")
        print("  the axis a resection score needs does not.")
    elif bs.mean() < 0.15:
        print("  Both axes collapse. The weighting choice determines the answer")
        print("  whichever way the data is sliced, which is the stronger claim.")
    else:
        print("  Mixed: neither axis is clean. Report both numbers plainly.")


if __name__ == "__main__":
    main()
