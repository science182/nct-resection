"""How often does a published-style conclusion change when the weighting does?

The rest of this repository establishes that weighting choice can determine a
degree-corrected map. That is a caution. This asks what the caution costs, by
taking the kinds of statement papers actually make and measuring how often each
one flips when only the edge weighting convention changes.

Five statement types, chosen because they appear routinely in structural
connectomics:

    top hub          "parcel X is the most central node"
    top-10 set       "these are the ten most damaging parcels to remove"
    rank order       the full ranking, compared by Spearman
    pairwise call    "removing A is worse than removing B"
    network claim    "network N is over-represented among high-risk parcels"

A flip rate near zero would mean the weighting debate is academic. A high flip
rate means published conclusions of that type are contingent on a preprocessing
choice that is usually reported in one clause of a methods section.

Everything runs from cached deletion sweeps, so no long recomputation.

    python3 run_fliprate.py
"""

import itertools
import os

import numpy as np
from scipy.stats import spearmanr

from run_confounds import residualize
from run_weightings import WEIGHTINGS

TOP_N = 10
PAIR_SAMPLES = 20000


def jaccard(a, b):
    sa, sb = set(a.tolist()), set(b.tolist())
    return len(sa & sb) / len(sa | sb)


def pair_flip_rate(x, y, rng, n_samples=PAIR_SAMPLES, margin=0.0):
    """Fraction of random parcel pairs whose ordering disagrees between maps.

    `margin` ignores pairs that are near-ties under the first map, since a
    coin-flip on two indistinguishable parcels is not a meaningful reversal.
    """
    n = len(x)
    i = rng.integers(0, n, n_samples)
    j = rng.integers(0, n, n_samples)
    keep = i != j
    i, j = i[keep], j[keep]

    dx = x[i] - x[j]
    if margin > 0:
        spread = np.percentile(np.abs(x - x.mean()), 75)
        keep = np.abs(dx) > margin * spread
        i, j, dx = i[keep], j[keep], dx[keep]
    dy = y[i] - y[j]
    return float(np.mean(np.sign(dx) != np.sign(dy))), len(i)


def conclusions(maps, labels, rng, asg=None):
    """Compare every pair of weightings on each statement type."""
    rows = []
    for a, b in itertools.combinations(labels, 2):
        x, y = maps[a], maps[b]
        top_x, top_y = np.argsort(-x)[:TOP_N], np.argsort(-y)[:TOP_N]
        flip_all, n_all = pair_flip_rate(x, y, rng)
        flip_clear, n_clear = pair_flip_rate(x, y, rng, margin=0.5)
        row = {
            "pair": f"{a} vs {b}",
            "hub_same": int(np.argmax(x) == np.argmax(y)),
            "top10_jaccard": jaccard(top_x, top_y),
            "spearman": float(spearmanr(x, y)[0]),
            "pair_flip": flip_all,
            "pair_flip_clear": flip_clear,
            "n_clear": n_clear,
        }
        if asg is not None:
            row["net_x"] = dominant_network(x, asg)
            row["net_y"] = dominant_network(y, asg)
            row["net_same"] = int(row["net_x"] == row["net_y"])
        rows.append(row)
    return rows


def dominant_network(score, asg, top=36):
    """Which network is most over-represented among the highest-risk parcels."""
    top_set = set(np.argsort(-score)[:top].tolist())
    best, best_ratio = None, 0.0
    for net in range(1, 11):
        members = set(np.flatnonzero((asg == net) | (asg == net + 10)).tolist())
        if not members:
            continue
        expected = top * len(members) / len(asg)
        observed = len(top_set & members)
        ratio = observed / expected if expected else 0.0
        if ratio > best_ratio:
            best, best_ratio = net, ratio
    return best


def summarize(name, rows, note=""):
    print(f"\n=== {name} ===")
    if note:
        print(f"  {note}")
    print(f"  {'weighting pair':<22} {'top hub':>8} {'top10':>7} {'rho':>7} "
          f"{'pair flips':>11} {'clear-cut':>10}")
    for r in rows:
        print(f"  {r['pair']:<22} {'same' if r['hub_same'] else 'DIFFER':>8} "
              f"{r['top10_jaccard']:>7.2f} {r['spearman']:>+7.2f} "
              f"{r['pair_flip']:>10.1%} {r['pair_flip_clear']:>10.1%}")

    hub = 100 * (1 - np.mean([r["hub_same"] for r in rows]))
    jac = np.mean([r["top10_jaccard"] for r in rows])
    flip = np.mean([r["pair_flip_clear"] for r in rows])
    print(f"  {'':<22} {'':>8} {'':>7} {'':>7} {'':>11} {'':>10}")
    print(f"  top hub changes in {hub:.0f}% of comparisons")
    print(f"  top-10 sets overlap {jac:.0%} on average")
    print(f"  clear-cut pairwise calls reverse {flip:.1%} of the time")
    if "net_same" in rows[0]:
        nets = 100 * (1 - np.mean([r["net_same"] for r in rows]))
        print(f"  dominant network changes in {nets:.0f}% of comparisons")


def load_maps(path, corrected):
    d = np.load(path)
    n = int(d["done"]) if "done" in d.files else d[f"{WEIGHTINGS[0]}_d_ac"].shape[0]
    out = {}
    for w in WEIGHTINGS:
        ac, st = d[f"{w}_d_ac"][:n], d[f"{w}_strength"][:n]
        out[w] = (residualize(ac, st, "quadratic").mean(0) if corrected
                  else ac.mean(0))
    return out


def main():
    rng = np.random.default_rng(0)
    try:
        from run_energy import load_network_assignment
        asg, _ = load_network_assignment()
    except Exception:
        asg = None

    print("\nHOW OFTEN DOES A CONCLUSION CHANGE WHEN ONLY THE WEIGHTING DOES?")
    print(f"  statements compared across {len(WEIGHTINGS)} conventions, "
          f"{len(list(itertools.combinations(WEIGHTINGS, 2)))} pairs each")

    hcp = "data/weighting_deletions.npz"
    if os.path.exists(hcp):
        raw = load_maps(hcp, corrected=False)
        summarize("HCP, raw deletion damage", conclusions(raw, WEIGHTINGS, rng, asg),
                  "the measure as usually reported")
        cor = load_maps(hcp, corrected=True)
        summarize("HCP, degree-corrected", conclusions(cor, WEIGHTINGS, rng, asg),
                  "after removing the degree component")
    else:
        print(f"\n  skipped HCP: {hcp} not present")

    lau = "data/lausanne_deletions.npz"
    if os.path.exists(lau):
        raw = load_maps(lau, corrected=False)
        summarize("Lausanne 219, raw deletion damage",
                  conclusions(raw, WEIGHTINGS, rng),
                  "independent dataset, milder weight distribution")
        cor = load_maps(lau, corrected=True)
        summarize("Lausanne 219, degree-corrected",
                  conclusions(cor, WEIGHTINGS, rng))
    else:
        print(f"\n  skipped Lausanne: {lau} not present")

    print("\n  Pairwise calls counted twice: over all random parcel pairs, and")
    print("  over 'clear-cut' pairs only, where the first map separates the two")
    print("  parcels by more than half its interquartile spread. The second is")
    print("  the honest number, since reversing a near-tie is not a real flip.")


if __name__ == "__main__":
    main()
