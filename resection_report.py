"""Convention-robust comparison of candidate resections.

A connectomic planning pipeline normally returns one number per candidate
resection: the damage it does to some network measure. That number depends on
choices nobody sees in the output. On a heavy-tailed connectome the
degree-corrected ranking can reorder entirely when the edge weighting changes,
and the published deletion designs use unthresholded streamline counts, which is
the heaviest-tailed case measured here.

The clinical question is rarely "how damaging is this parcel" in the abstract.
It is "is corridor A safer than corridor B". That comparison can be robust even
when the underlying numbers are unstable, and it can be a coin flip even when
each number looks precise. This reports which.

For every pair of candidates it asks how often A beats B across the analytic
choices, and grades the comparison:

    ROBUST        the same winner under every choice
    LIKELY        one winner in at least four fifths of choices
    COIN FLIP     the answer is set by the analysis, not the anatomy

The choices varied are the edge weighting convention (raw, log, row-normalized,
binarized) crossed with the damage measure (average controllability, global
efficiency, node strength). Strength is included deliberately: if it agrees with
everything else, the sophisticated measures are not earning their place.

    python3 resection_report.py --input conn.mat --candidates 12,45,88
    python3 resection_report.py --input conn.mat --grow 5 --seeds 12,45,88
"""

import argparse
import itertools

import numpy as np

from audit import clean, load_any, tail_ratio
from controllability import average_controllability, spectral_scale
from lesion import global_efficiency
from run_weightings import WEIGHTINGS, reweight

MEASURES = ("controllability", "efficiency", "strength")


def damage(A, resected, measure):
    """Damage done by removing `resected`, scored the requested way.

    Higher means worse. Controllability and efficiency are scored over the
    surviving parcels so the comparison is not confounded by network size.
    """
    n = A.shape[0]
    keep = np.setdiff1d(np.arange(n), np.asarray(resected, dtype=int))
    if keep.size < 2:
        return float("inf")

    if measure == "strength":
        return float(A.sum(axis=1)[np.asarray(resected, dtype=int)].sum())
    if measure == "efficiency":
        return float(global_efficiency(A) - global_efficiency(A, resected))

    ref = spectral_scale(A)
    intact = average_controllability(A, ref_scale=ref)[keep].mean()
    lesioned = average_controllability(
        A[np.ix_(keep, keep)], ref_scale=ref).mean()
    return float(intact - lesioned)


def score_candidates(A, candidates, weightings=WEIGHTINGS, measures=MEASURES):
    """damage[(weighting, measure)] -> array over candidates."""
    scores = {}
    for w in weightings:
        Aw = reweight(A, w)
        for m in measures:
            scores[(w, m)] = np.array([damage(Aw, c, m) for c in candidates])
    return scores


def grade_pairs(scores, candidates, labels=None):
    """For each pair, how consistently does one candidate beat the other."""
    labels = labels or [str(sorted(c)) for c in candidates]
    n = len(candidates)
    rows = []
    for i, j in itertools.combinations(range(n), 2):
        # Ties must not count as wins. Two candidates that produce identical
        # damage under every choice would otherwise be graded ROBUST in favour
        # of whichever came second, which is how identical corridors were
        # reported as a reliable difference.
        wins_i = sum(1 for s in scores.values() if s[i] < s[j])
        ties = sum(1 for s in scores.values() if s[i] == s[j])
        total = len(scores)
        frac = (wins_i + 0.5 * ties) / total
        safer, share = (i, frac) if frac >= 0.5 else (j, 1 - frac)
        if ties == total:
            grade = "IDENTICAL"
        elif share == 1.0:
            grade = "ROBUST"
        elif share >= 0.8:
            grade = "LIKELY"
        else:
            grade = "COIN FLIP"
        rows.append({
            "a": i, "b": j, "safer": safer, "share": share, "grade": grade,
            "label_a": labels[i], "label_b": labels[j],
            "label_safer": labels[safer],
        })
    return rows


def consensus_rank(scores, n_candidates):
    """Mean rank across every analytic choice, plus how much that rank moves."""
    ranks = np.zeros((len(scores), n_candidates))
    for k, s in enumerate(scores.values()):
        order = np.argsort(s)
        r = np.empty(n_candidates)
        r[order] = np.arange(1, n_candidates + 1)
        ranks[k] = r
    return ranks.mean(0), ranks.std(0), ranks.min(0), ranks.max(0)


def parse_candidates(text):
    """'1,2,3' -> one candidate; '1,2|7,8' -> two multi-parcel candidates."""
    out = []
    for chunk in text.split("|"):
        out.append([int(x) for x in chunk.split(",") if x.strip() != ""])
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--log10", action="store_true")
    ap.add_argument("--subject", type=int, default=0,
                    help="which matrix to use if the file holds several")
    ap.add_argument("--candidates",
                    help="parcel indices; comma within a candidate, | between "
                         "candidates, e.g. 12,13|45,46")
    ap.add_argument("--seeds", help="with --grow, seed parcels to grow from")
    ap.add_argument("--grow", type=int, default=0,
                    help="grow each seed into a contiguous resection of this size")
    ap.add_argument("--max-nodes", type=int, default=400)
    args = ap.parse_args()

    A = clean(load_any(args.input, log10=args.log10,
                       max_nodes=args.max_nodes)[args.subject])
    n = A.shape[0]

    if args.grow:
        if not args.seeds:
            ap.error("--grow needs --seeds")
        from lesion import grow_resection
        adjacency = None
        if n == 360:
            try:
                from annot import build_hcp_mmp1_adjacency
                adjacency = build_hcp_mmp1_adjacency()
            except Exception:
                adjacency = None
        seeds = [int(x) for x in args.seeds.split(",")]
        candidates = [grow_resection(A, s, args.grow, adjacency=adjacency)
                      for s in seeds]
        if adjacency is None:
            print("note: no parcel adjacency available, corridors grown along"
                  " connectivity weight instead\n")
    elif args.candidates:
        candidates = parse_candidates(args.candidates)
    else:
        ap.error("give --candidates or --seeds with --grow")

    for c in candidates:
        if max(c) >= n or min(c) < 0:
            ap.error(f"parcel index out of range for a {n}-parcel connectome: {c}")

    # Growing from different seeds can converge on the same corridor, which
    # then gets compared against itself.
    seen, unique, dropped = {}, [], []
    for c in candidates:
        key = tuple(sorted(c))
        if key in seen:
            dropped.append(key)
        else:
            seen[key] = True
            unique.append(c)
    if dropped:
        print(f"note: {len(dropped)} duplicate candidate(s) removed; "
              f"different seeds grew into the same corridor")
        for key in dropped:
            print(f"      {','.join(str(i) for i in key)}")
        print()
    candidates = unique
    if len(candidates) < 2:
        ap.error("need at least two distinct candidates to compare")

    labels = [",".join(str(i) for i in sorted(c)) for c in candidates]
    scores = score_candidates(A, candidates)
    n_choices = len(scores)

    print(f"\nRESECTION COMPARISON: {args.input}")
    print(f"  {n} parcels, {len(candidates)} candidates, "
          f"{n_choices} analytic choices "
          f"({len(WEIGHTINGS)} weightings x {len(MEASURES)} measures)")
    print(f"  tail ratio {tail_ratio(A):.0f}"
          f"{'   HEAVY, expect convention sensitivity' if tail_ratio(A) > 200 else ''}\n")

    mean_r, sd_r, min_r, max_r = consensus_rank(scores, len(candidates))
    print("CONSENSUS RANKING  (1 = least damaging)")
    print(f"  {'candidate':<26} {'mean rank':>10} {'range':>10} {'stability':>11}")
    for i in np.argsort(mean_r):
        span = f"{int(min_r[i])}-{int(max_r[i])}"
        stable = "firm" if sd_r[i] < 0.5 else "soft" if sd_r[i] < 1.5 else "unstable"
        print(f"  {labels[i]:<26} {mean_r[i]:>10.2f} {span:>10} {stable:>11}")

    print("\nPAIRWISE DECISIONS")
    rows = grade_pairs(scores, candidates, labels)
    for r in sorted(rows, key=lambda r: -r["share"]):
        if r["grade"] == "COIN FLIP":
            print(f"  [{r['grade']:^9}] {r['label_a']} vs {r['label_b']}: "
                  f"no reliable winner ({r['share']:.0%} agreement)")
        else:
            print(f"  [{r['grade']:^9}] {r['label_safer']} is safer "
                  f"({r['share']:.0%} of choices)")

    robust = sum(1 for r in rows if r["grade"] == "ROBUST")
    flips = sum(1 for r in rows if r["grade"] == "COIN FLIP")
    print(f"\n  {robust} of {len(rows)} comparisons robust, {flips} coin flips")

    agree = np.mean([
        np.corrcoef(scores[(w, "strength")], scores[(w, m)])[0, 1]
        for w in WEIGHTINGS for m in ("controllability", "efficiency")
        if np.std(scores[(w, m)]) > 0 and np.std(scores[(w, "strength")]) > 0])
    print(f"  agreement with plain node strength: {agree:+.2f}")
    if agree > 0.9:
        print("  The network measures are tracking degree here. Reporting the")
        print("  total connectivity removed would say the same thing.")


if __name__ == "__main__":
    main()
