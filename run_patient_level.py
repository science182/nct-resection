"""How stable is a conclusion for ONE patient, not for a cohort average?

Every stability figure elsewhere in this repository is computed on a map averaged
over subjects. No surgeon sees a cohort average. They see one brain, once, and
make a decision from it.

Averaging suppresses noise, so a group map is stable by construction in a way a
single subject's is not. If convention-sensitivity is already substantial on the
average, the per-patient figure is the one that matters for a planning product,
and it should be worse. This measures how much worse.

Three quantities, all per subject and then summarized across subjects:

    within-patient stability   for one subject, correlate that subject's own map
                               under one weighting against their own map under
                               another. This is what a clinician is exposed to.

    group stability            the same comparison on the cohort-mean map, for
                               reference. Reported elsewhere in this repo.

    between-patient spread     under a fixed weighting, how much do subjects
                               disagree with each other? This bounds how much of
                               the per-patient instability is weighting rather
                               than ordinary individual variation.

The third is the control that makes the first interpretable. If subjects disagree
with each other as much as a subject disagrees with itself across weightings,
then weighting is not the dominant problem and individual variation is.

    python3 run_patient_level.py
"""

import itertools
import os

import numpy as np
from scipy.stats import spearmanr

from run_confounds import residualize
from run_fliprate import TOP_N, jaccard
from run_weightings import WEIGHTINGS

DEL = "data/weighting_deletions.npz"
MONOTONE = ("raw", "log", "rownorm")


def load(path, corrected):
    d = np.load(path)
    n = int(d["done"]) if "done" in d.files else d["raw_d_ac"].shape[0]
    out = {}
    for w in WEIGHTINGS:
        ac, st = d[f"{w}_d_ac"][:n], d[f"{w}_strength"][:n]
        out[w] = residualize(ac, st, "quadratic") if corrected else ac
    return out, n


def within_patient(maps, conventions, n):
    """Per subject, agreement between that subject's own maps."""
    rho, jac = [], []
    for s in range(n):
        rs, js = [], []
        for a, b in itertools.combinations(conventions, 2):
            x, y = maps[a][s], maps[b][s]
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            rs.append(spearmanr(x, y)[0])
            js.append(jaccard(np.argsort(-x)[:TOP_N], np.argsort(-y)[:TOP_N]))
        if rs:
            rho.append(np.mean(rs))
            jac.append(np.mean(js))
    return np.array(rho), np.array(jac)


def group_level(maps, conventions):
    """The same comparison on cohort-mean maps."""
    rho, jac = [], []
    for a, b in itertools.combinations(conventions, 2):
        x, y = maps[a].mean(0), maps[b].mean(0)
        rho.append(spearmanr(x, y)[0])
        jac.append(jaccard(np.argsort(-x)[:TOP_N], np.argsort(-y)[:TOP_N]))
    return float(np.mean(rho)), float(np.mean(jac))


def between_patient(maps, conventions, n, rng):
    """Under a fixed weighting, how much do different subjects disagree?"""
    rho, jac = [], []
    pairs = [(i, j) for i, j in itertools.combinations(range(n), 2)]
    if len(pairs) > 200:
        idx = rng.choice(len(pairs), 200, replace=False)
        pairs = [pairs[k] for k in idx]
    for w in conventions:
        for i, j in pairs:
            x, y = maps[w][i], maps[w][j]
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            rho.append(spearmanr(x, y)[0])
            jac.append(jaccard(np.argsort(-x)[:TOP_N], np.argsort(-y)[:TOP_N]))
    return float(np.mean(rho)), float(np.mean(jac))


def report(label, corrected, conventions):
    maps, n = load(DEL, corrected)
    rng = np.random.default_rng(0)

    w_rho, w_jac = within_patient(maps, conventions, n)
    g_rho, g_jac = group_level(maps, conventions)
    b_rho, b_jac = between_patient(maps, conventions, n, rng)

    print(f"\n=== {label} ===")
    print(f"  {n} subjects, conventions: {', '.join(conventions)}")
    print(f"  {'comparison':<40} {'rho':>8} {'top-10 overlap':>15}")
    print(f"  {'one patient, across weightings':<40} {w_rho.mean():>+8.3f} "
          f"{w_jac.mean():>15.0%}")
    print(f"  {'cohort mean, across weightings':<40} {g_rho:>+8.3f} "
          f"{g_jac:>15.0%}")
    print(f"  {'two patients, same weighting':<40} {b_rho:>+8.3f} "
          f"{b_jac:>15.0%}")

    print(f"  per-patient spread: rho ranges {w_rho.min():+.3f} to "
          f"{w_rho.max():+.3f}, sd {w_rho.std():.3f}")
    worse = g_rho - w_rho.mean()
    print(f"  averaging inflates apparent stability by {worse:+.3f} in rho, "
          f"{g_jac - w_jac.mean():+.0%} in overlap")
    if w_rho.mean() < b_rho:
        print("  A patient agrees with another patient MORE than with their own")
        print("  map under a different weighting. The analytic choice moves the")
        print("  result further than the biology does.")
    else:
        print("  Individual variation exceeds the weighting effect at the")
        print("  patient level, so weighting is not the dominant term here.")
    return w_rho, w_jac


def main():
    if not os.path.exists(DEL):
        print(f"needs {DEL}; run run_weightings.py first")
        return

    print("\nPER-PATIENT STABILITY, WHICH IS WHAT A PLANNING PRODUCT FACES")
    print("  Group-average figures elsewhere in this repo suppress subject")
    print("  noise by averaging. A clinician sees one brain.")

    report("raw deletion damage, monotone conventions", False, MONOTONE)
    report("degree-corrected, monotone conventions", True, MONOTONE)
    report("raw deletion damage, all four conventions", False, WEIGHTINGS)


if __name__ == "__main__":
    main()
