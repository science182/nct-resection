"""Is any resection measure stable across edge-weighting conventions?

Four weightings derived from the same streamline counts, so the tractography is
identical and only the convention changes. Two things to find out.

Mechanism: Fpt is streamline counts normalized per seed parcel. If `rownorm`
reproduces the language peak that Fpt gave, then per-seed normalization alone
explains why the two datasets disagreed, and the disagreement is fully
understood rather than mysterious.

Stability: a score whose spatial map changes when you change convention cannot
support a clinical risk score, whatever its p-values say. The comparison that
matters is whether global efficiency, the measure this project set out to beat,
is any steadier than controllability under the same test.
"""

import numpy as np
from scipy.stats import spearmanr

from run_confounds import K, residualize
from run_energy import load_network_assignment
from run_weightings import OUT, WEIGHTINGS
from spin_test import parcel_centroids, spin_permutations

NETWORKS = {"Language": 5, "Visual": 1}


def enrichment(score, asg, perms, net):
    members = set(np.flatnonzero((asg == net) | (asg == net + 10)).tolist())
    obs = len(set(np.argsort(-score)[:K].tolist()) & members)
    spun = np.array([len(set(np.argsort(-score[p])[:K].tolist()) & members)
                     for p in perms])
    return obs, (spun >= obs).mean(), K * len(members) / 360


def main():
    d = np.load(OUT)
    n_done = int(d["done"])
    asg, _ = load_network_assignment()
    print(f"{n_done} subjects completed\n")

    print("building spin permutations...")
    perms = spin_permutations(parcel_centroids("lh", "L"),
                              parcel_centroids("rh", "R"), n_spins=1000)

    maps = {}
    for w in WEIGHTINGS:
        for meas in ("d_ac", "d_ge"):
            prof = residualize(d[f"{w}_{meas}"][:n_done],
                               d[f"{w}_strength"][:n_done], "quadratic")
            maps[(w, meas)] = prof.mean(0)

    print("\n=== 1. what does each weighting flag? (spin-tested) ===")
    print(f"  {'weighting':<10} {'measure':<7} "
          f"{'Language':>18} {'Visual':>18}")
    for w in WEIGHTINGS:
        for meas in ("d_ac", "d_ge"):
            cells = []
            for name, net in NETWORKS.items():
                obs, p, exp = enrichment(maps[(w, meas)], asg, perms, net)
                cells.append(f"{obs:>2}/{K} p={p:.3f}")
            print(f"  {w:<10} {meas:<7} {cells[0]:>18} {cells[1]:>18}")
    print(f"  (chance: Language 2.3, Visual 6.0 of {K})")

    print("\n=== 2. do the maps agree across weightings? ===")
    for meas, label in (("d_ac", "controllability"), ("d_ge", "global efficiency")):
        print(f"  {label} residual, pairwise Spearman:")
        vals = []
        for i, a in enumerate(WEIGHTINGS):
            row = []
            for b in WEIGHTINGS:
                r = spearmanr(maps[(a, meas)], maps[(b, meas)])[0]
                row.append(f"{r:+.2f}")
                if a < b:
                    vals.append(r)
            print(f"    {a:<9} " + "  ".join(row))
        print(f"    mean off-diagonal agreement: {np.mean(vals):+.3f}\n")

    print("=== 3. verdict ===")
    ac_stab = np.mean([spearmanr(maps[(a, 'd_ac')], maps[(b, 'd_ac')])[0]
                       for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]])
    ge_stab = np.mean([spearmanr(maps[(a, 'd_ge')], maps[(b, 'd_ge')])[0]
                       for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]])
    print(f"  controllability map stability across conventions: {ac_stab:+.3f}")
    print(f"  global efficiency map stability across conventions: {ge_stab:+.3f}")
    steadier = "global efficiency" if ge_stab > ac_stab else "controllability"
    print(f"  steadier measure: {steadier}")
    if max(ac_stab, ge_stab) < 0.5:
        print("  neither is stable enough to carry a clinical score as it stands")


if __name__ == "__main__":
    main()
