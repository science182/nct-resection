"""Is there a resection score sensitive to position rather than amount?

Every deletion measure tested so far inherits "how much connectivity did you
remove", because cutting a parcel removes exactly its own strength. Average
controllability tracks strength at +0.892 within subject for that reason.

So define the score as what is left after removing that dependence, and then
ask the only question that matters: is the leftover a real, reproducible
spatial pattern, or is it noise?

Four tests:

1. Split-half reproducibility. Split the 24 subjects, average the residual
   profile in each half, correlate. This is the decisive one. A reproducible
   spatial pattern means position carries real signal. Near zero means the
   residual is noise and there is nothing here.

2. Does the controllability residual differ from the efficiency residual? If
   they agree, controllability adds nothing that global efficiency does not
   already provide once both are degree-corrected.

3. Strength-matched pairs. Among parcels with nearly identical strength,
   does damage still differ systematically and consistently across subjects?
   This tests position sensitivity without relying on regression at all.

4. Anatomical organization. Does the residual align with network membership,
   or is it spatially arbitrary?
"""

import numpy as np
from scipy.stats import spearmanr

from run_energy import load_network_assignment

DEL = "data/subject_deletions.npz"


def residualize_within_subject(damage, strength):
    """Remove each subject's own strength dependence, parcel-wise across the brain.

    Done per subject rather than per parcel: the question is which parcels are
    riskier than their connectivity alone predicts, inside that person's brain.
    Rank-based, since the strength-damage relationship need not be linear.
    """
    out = np.zeros_like(damage)
    for s in range(damage.shape[0]):
        x = np.argsort(np.argsort(strength[s])).astype(float)
        y = np.argsort(np.argsort(damage[s])).astype(float)
        design = np.column_stack([np.ones_like(x), x, x ** 2])
        out[s] = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return out


def split_half(profiles, n_rep=200, seed=0):
    """Correlate mean residual profiles between random halves of subjects."""
    rng = np.random.default_rng(seed)
    n = profiles.shape[0]
    rs = []
    for _ in range(n_rep):
        perm = rng.permutation(n)
        a, b = perm[: n // 2], perm[n // 2:]
        rs.append(spearmanr(profiles[a].mean(0), profiles[b].mean(0))[0])
    return np.array(rs)


def main():
    d = np.load(DEL)
    d_ac, d_ge, st = d["d_ac"], d["d_ge"], d["strength"]
    n_subj, n_parcel = d_ac.shape
    print(f"{n_subj} subjects, {n_parcel} parcels\n")

    res_ac = residualize_within_subject(d_ac, st)
    res_ge = residualize_within_subject(d_ge, st)

    print("=== 0. did residualizing actually remove strength? ===")
    for name, r in (("AC residual", res_ac), ("GE residual", res_ge)):
        rho = np.mean([spearmanr(r[s], st[s])[0] for s in range(n_subj)])
        print(f"  {name} vs strength, within subject: {rho:+.4f}")

    print("\n=== 1. split-half reproducibility (the decisive test) ===")
    for name, r in (("AC residual", res_ac), ("GE residual", res_ge),
                    ("raw d_ac", d_ac), ("strength", st)):
        rs = split_half(r)
        print(f"  {name:<14} r = {rs.mean():+.3f}  (sd {rs.std():.3f}, "
              f"2.5-97.5 pct {np.percentile(rs, 2.5):+.3f} to "
              f"{np.percentile(rs, 97.5):+.3f})")
    print("  raw measures set the ceiling; the residual rows are the question")

    print("\n=== 2. is the AC residual different from the GE residual? ===")
    per_subj = np.mean([spearmanr(res_ac[s], res_ge[s])[0] for s in range(n_subj)])
    group = spearmanr(res_ac.mean(0), res_ge.mean(0))[0]
    print(f"  within subject: {per_subj:+.3f}")
    print(f"  group mean profiles: {group:+.3f}")
    print(f"  -> {'redundant with GE' if abs(group) > 0.8 else 'carries distinct information'}")

    print("\n=== 3. strength-matched pairs ===")
    mean_st = st.mean(0)
    mean_res = res_ac.mean(0)
    order = np.argsort(mean_st)
    diffs, gaps = [], []
    for k in range(0, len(order) - 1, 2):
        i, j = order[k], order[k + 1]
        rel_gap = abs(mean_st[i] - mean_st[j]) / mean_st[[i, j]].mean()
        if rel_gap < 0.02:  # within 2 percent strength
            diffs.append(abs(mean_res[i] - mean_res[j]))
            gaps.append(rel_gap)
    diffs = np.array(diffs)
    spread = np.abs(mean_res[:, None] - mean_res[None, :])
    spread = spread[~np.eye(n_parcel, dtype=bool)]
    print(f"  {len(diffs)} parcel pairs matched within 2% on strength")
    print(f"  median |residual difference| within matched pairs: {np.median(diffs):.1f}")
    print(f"  median |residual difference| between random parcels: {np.median(spread):.1f}")
    print(f"  ratio: {np.median(diffs) / np.median(spread):.2f}  "
          f"(near 1 means strength matching does not reduce the difference, "
          f"so position drives it)")

    print("\n=== 4. does the residual respect network anatomy? ===")
    asg, labels = load_network_assignment()
    within, between = [], []
    for i in range(n_parcel):
        for j in range(i + 1, n_parcel):
            dv = abs(mean_res[i] - mean_res[j])
            (within if asg[i] == asg[j] else between).append(dv)
    within, between = np.mean(within), np.mean(between)
    print(f"  mean |difference| within the same network:  {within:.1f}")
    print(f"  mean |difference| across different networks: {between:.1f}")
    print(f"  ratio: {within / between:.3f}  "
          f"({'organized by network' if within / between < 0.95 else 'no network organization'})")


if __name__ == "__main__":
    main()
