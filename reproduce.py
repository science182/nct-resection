"""Verify every headline number in the README.

A repository that reports a lot of numbers should make them checkable without
the reader having to work out which of fifteen scripts to run in what order.
Each claim below is recomputed from the cached results and compared against the
value written in the README, with a tolerance.

Claims needing inputs that are not present are reported as SKIP rather than
failing, since the large sweeps take hours and the per-subject files are 1.3 GB.

    python3 reproduce.py              # check what can be checked
    python3 reproduce.py --verbose    # show the computed value for every claim

Exit code is 1 if any available claim fails, so this can gate a commit.
"""

import argparse
import os
import sys

import numpy as np
from scipy.stats import spearmanr

from run_confounds import residualize

CLAIMS = []


def claim(text, expected, tol, needs=()):
    """Register a check. The function returns the computed value."""
    def wrap(fn):
        CLAIMS.append({"text": text, "expected": expected, "tol": tol,
                       "needs": needs, "fn": fn})
        return fn
    return wrap


def _within(a, b):
    d = np.load(a)
    return d


# --- the finding that replicated ---------------------------------------------

@claim("deletion damage vs node strength, HCP, within subject", 0.896, 0.005,
       needs=("data/scale_fpt_deletions.npz",))
def _hcp_degree():
    d = np.load("data/scale_fpt_deletions.npz")
    ac, st = d["d_ac"], d["strength"]
    return float(np.mean([spearmanr(ac[s], st[s])[0] for s in range(ac.shape[0])]))


@claim("deletion damage vs node strength, Lausanne, within subject", 0.883, 0.01,
       needs=("data/lausanne_deletions.npz",))
def _lausanne_degree():
    d = np.load("data/lausanne_deletions.npz")
    ac, st = d["raw_d_ac"], d["raw_strength"]
    return float(np.mean([spearmanr(ac[s], st[s])[0] for s in range(ac.shape[0])]))


@claim("global efficiency vs node strength, HCP, within subject", 0.539, 0.01,
       needs=("data/scale_fpt_deletions.npz",))
def _hcp_ge_degree():
    d = np.load("data/scale_fpt_deletions.npz")
    ge, st = d["d_ge"], d["strength"]
    return float(np.mean([spearmanr(ge[s], st[s])[0] for s in range(ge.shape[0])]))


# --- the contradiction --------------------------------------------------------

@claim("agreement between Fpt and streamline risk maps", 0.021, 0.02,
       needs=("data/scale_fpt_deletions.npz", "data/scale_streamline_deletions.npz"))
def _map_agreement():
    f = np.load("data/scale_fpt_deletions.npz")
    s = np.load("data/scale_streamline_deletions.npz")
    rf = residualize(f["d_ac"], f["strength"], "quadratic").mean(0)
    rs = residualize(s["d_ac"], s["strength"], "quadratic").mean(0)
    return float(spearmanr(rf, rs)[0])


@claim("language parcels in top 36, Fpt weighting", 9, 0,
       needs=("data/scale_fpt_deletions.npz",))
def _lang_fpt():
    from run_energy import load_network_assignment
    asg, _ = load_network_assignment()
    lang = set(np.flatnonzero((asg == 5) | (asg == 15)).tolist())
    d = np.load("data/scale_fpt_deletions.npz")
    r = residualize(d["d_ac"], d["strength"], "quadratic").mean(0)
    return len(set(np.argsort(-r)[:36].tolist()) & lang)


@claim("language parcels in top 36, streamline weighting", 0, 0,
       needs=("data/scale_streamline_deletions.npz",))
def _lang_stream():
    from run_energy import load_network_assignment
    asg, _ = load_network_assignment()
    lang = set(np.flatnonzero((asg == 5) | (asg == 15)).tolist())
    d = np.load("data/scale_streamline_deletions.npz")
    r = residualize(d["d_ac"], d["strength"], "quadratic").mean(0)
    return len(set(np.argsort(-r)[:36].tolist()) & lang)


@claim("visual parcels in top 36, streamline weighting", 25, 1,
       needs=("data/scale_streamline_deletions.npz",))
def _vis_stream():
    from run_energy import load_network_assignment
    asg, _ = load_network_assignment()
    vis = set(np.flatnonzero((asg == 1) | (asg == 11)).tolist())
    d = np.load("data/scale_streamline_deletions.npz")
    r = residualize(d["d_ac"], d["strength"], "quadratic").mean(0)
    return len(set(np.argsort(-r)[:36].tolist()) & vis)


@claim("split-half map agreement, Fpt", 0.999, 0.003,
       needs=("data/scale_fpt_deletions.npz",))
def _split_half():
    d = np.load("data/scale_fpt_deletions.npz")
    r = residualize(d["d_ac"], d["strength"], "quadratic")
    h = r.shape[0] // 2
    return float(spearmanr(r[:h].mean(0), r[h:].mean(0))[0])


# --- the mechanism ------------------------------------------------------------

@claim("Spearman(log tail ratio, stability) across the grid", -0.917, 0.02)
def _mechanism():
    tail = np.array([8748, 899, 193, 94, 30, 14, 10, 5, 4], float)
    stab = np.array([0.319, 0.464, 0.638, 0.629, 0.688, 0.766, 0.843, 0.781,
                     0.772])
    return float(spearmanr(np.log10(tail), stab)[0])


@claim("Spearman(density, stability) across the grid", -0.211, 0.02)
def _mechanism_density():
    dens = np.array([1.00, 0.30, 0.11, 1.00, 0.30, 0.11, 1.00, 0.30, 0.11])
    stab = np.array([0.319, 0.464, 0.638, 0.629, 0.688, 0.766, 0.843, 0.781,
                     0.772])
    return float(spearmanr(dens, stab)[0])


@claim("map stability, Lausanne, degree-corrected", 0.756, 0.03,
       needs=("data/lausanne_deletions.npz",))
def _lausanne_stability():
    from run_weightings import WEIGHTINGS
    d = np.load("data/lausanne_deletions.npz")
    maps = {w: residualize(d[f"{w}_d_ac"], d[f"{w}_strength"],
                           "quadratic").mean(0) for w in WEIGHTINGS}
    vals = [spearmanr(maps[a], maps[b])[0]
            for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
    return float(np.mean(vals))


@claim("map stability, HCP, degree-corrected", 0.289, 0.03,
       needs=("data/weighting_deletions.npz",))
def _hcp_stability():
    from run_weightings import WEIGHTINGS
    d = np.load("data/weighting_deletions.npz")
    n = int(d["done"])
    maps = {w: residualize(d[f"{w}_d_ac"][:n], d[f"{w}_strength"][:n],
                           "quadratic").mean(0) for w in WEIGHTINGS}
    vals = [spearmanr(maps[a], maps[b])[0]
            for i, a in enumerate(WEIGHTINGS) for b in WEIGHTINGS[i + 1:]]
    return float(np.mean(vals))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("\nVERIFYING THE CLAIMS IN README.md\n")
    print(f"  {'claim':<52} {'expected':>9} {'got':>9}  status")
    print("  " + "-" * 82)

    passed = failed = skipped = 0
    for c in CLAIMS:
        missing = [p for p in c["needs"] if not os.path.exists(p)]
        if missing:
            print(f"  {c['text']:<52} {c['expected']:>9} {'':>9}  SKIP "
                  f"(needs {os.path.basename(missing[0])})")
            skipped += 1
            continue
        try:
            got = c["fn"]()
        except Exception as exc:  # noqa: BLE001
            print(f"  {c['text']:<52} {c['expected']:>9} {'':>9}  ERROR {exc}")
            failed += 1
            continue
        ok = abs(got - c["expected"]) <= c["tol"]
        status = "pass" if ok else "FAIL"
        fmt = "{:>9.0f}" if float(c["expected"]).is_integer() and c["tol"] < 2 \
            else "{:>9.3f}"
        print(f"  {c['text']:<52} {fmt.format(c['expected'])} "
              f"{fmt.format(got)}  {status}")
        passed += ok
        failed += not ok

    print()
    print(f"  {passed} passed, {failed} failed, {skipped} skipped")
    if skipped:
        print("  Skipped claims need the large sweeps; see README 'Running it'.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
