"""Rank simulated resections by each metric and find where the metrics disagree.

The disagreement set is the actual scientific output. If controllability ranks
the same resections that global efficiency already flags, controllability adds
nothing here and the project should stop. If it flags a distinct set, those
resections are the ones worth characterizing.
"""

import numpy as np
from scipy.stats import spearmanr

from controllability import (
    average_controllability,
    delta_controllability,
    modal_controllability,
)
from lesion import delta_global_efficiency, grow_resection, pagerank_scores


def sweep_resections(A, size, seeds=None, adjacency=None):
    """Simulate one contiguous resection per seed and score it every way.

    Returns a list of dicts, one per resection.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    seeds = range(n) if seeds is None else seeds
    pr = pagerank_scores(A)

    rows = []
    for seed in seeds:
        resected = grow_resection(A, seed, size, adjacency=adjacency)
        if len(resected) < size:
            continue  # could not grow a contiguous resection of the requested size
        try:
            ac = delta_controllability(A, resected, metric=average_controllability)
            mc = delta_controllability(A, resected, metric=modal_controllability)
            dge = delta_global_efficiency(A, resected)
        except ValueError:
            continue

        rows.append(
            {
                "seed": int(seed),
                "resected": resected,
                "d_ave_control": ac["drop"],
                "d_mod_control": mc["drop"],
                "d_global_eff": dge,
                "seed_pagerank": float(pr[seed]),
            }
        )
    return rows


def rank_desc(values):
    """Rank so that 1 = most damaging."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(-values)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks


def disagreement(rows, key_a="d_ave_control", key_b="d_global_eff", top_frac=0.1):
    """Resections that one metric puts in its worst tier and the other does not."""
    a = rank_desc([r[key_a] for r in rows])
    b = rank_desc([r[key_b] for r in rows])
    cutoff = max(1, int(round(len(rows) * top_frac)))

    only_a = [rows[i] for i in range(len(rows)) if a[i] <= cutoff and b[i] > cutoff]
    only_b = [rows[i] for i in range(len(rows)) if b[i] <= cutoff and a[i] > cutoff]

    rho, p = spearmanr([r[key_a] for r in rows], [r[key_b] for r in rows])
    return {
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "cutoff_n": cutoff,
        f"top_{key_a}_only": only_a,
        f"top_{key_b}_only": only_b,
    }


def summarize(rows, top_frac=0.1):
    lines = []
    d = disagreement(rows, top_frac=top_frac)

    lines.append(f"resections scored: {len(rows)}")
    lines.append(
        f"Spearman(delta average controllability, delta global efficiency) = "
        f"{d['spearman_rho']:.3f}  (p = {d['spearman_p']:.2g})"
    )
    lines.append("")
    lines.append(
        f"Top {int(top_frac * 100)}% most damaging by each metric, "
        f"n = {d['cutoff_n']} per tier"
    )

    only_ac = d["top_d_ave_control_only"]
    only_ge = d["top_d_global_eff_only"]
    lines.append(
        f"  flagged by controllability but NOT by global efficiency: "
        f"{[r['seed'] for r in only_ac]}"
    )
    lines.append(
        f"  flagged by global efficiency but NOT by controllability: "
        f"{[r['seed'] for r in only_ge]}"
    )
    lines.append("")

    if not only_ac and not only_ge:
        lines.append(
            "  The metrics agree completely at this threshold. On real data that "
            "result would mean controllability is redundant with GE."
        )
    else:
        lines.append(
            "  These are the resections worth characterizing: they are the only "
            "place controllability could add clinical information."
        )

    # Does PageRank predict the worst resection? This mirrors the 35-75% figure
    # reported in Lin et al. 2024 and is the residual the proposal targets.
    for key, label in [
        ("d_global_eff", "global efficiency"),
        ("d_ave_control", "average controllability"),
    ]:
        worst = max(rows, key=lambda r: r[key])
        by_pr = sorted(rows, key=lambda r: -r["seed_pagerank"])
        hit_top1 = by_pr[0]["seed"] == worst["seed"]
        hit_top3 = worst["seed"] in [r["seed"] for r in by_pr[:3]]
        lines.append(
            f"  highest-PageRank seed is the worst resection by {label}: "
            f"top-1 {hit_top1}, within top-3 {hit_top3}"
        )

    return "\n".join(lines)
