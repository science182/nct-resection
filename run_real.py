"""Single-parcel deletion sweep on the real group-average connectome.

Single-parcel rather than contiguous multi-parcel, deliberately: the Rosen &
Halgren matrix is fully dense, so connectivity weight is useless as a proxy for
anatomical adjacency and any "contiguous" resection it grows would be
anatomically meaningless. Multi-parcel work waits for real parcel adjacency
derived from the HCP-MMP1 annotation files.

Thresholding is swept because a dense probabilistic connectome spans six orders
of magnitude in edge weight, and every ranking below depends on where the cut
is made. If the disagreement between controllability and global efficiency only
appears at one threshold, that is a finding about the threshold, not the brain.
"""

import numpy as np
from scipy.stats import spearmanr

from controllability import (
    average_controllability,
    modal_controllability,
    spectral_scale,
)
from data import load_rosen_halgren
from lesion import global_efficiency, pagerank_scores

CSV = "data/averageConnectivity_Fpt.csv"


def sweep_single_parcel(A):
    """Delete each parcel in turn, scoring every metric. Intact quantities are
    computed once rather than per deletion."""
    n = A.shape[0]
    ref = spectral_scale(A)
    intact_ac = average_controllability(A, ref_scale=ref)
    intact_mc = modal_controllability(A, ref_scale=ref)
    baseline_ge = global_efficiency(A)
    pr = pagerank_scores(A)

    d_ac = np.zeros(n)
    d_mc = np.zeros(n)
    d_ge = np.zeros(n)

    for i in range(n):
        retained = np.setdiff1d(np.arange(n), [i])
        sub = A[np.ix_(retained, retained)]
        d_ac[i] = intact_ac[retained].mean() - average_controllability(
            sub, ref_scale=ref).mean()
        d_mc[i] = intact_mc[retained].mean() - modal_controllability(
            sub, ref_scale=ref).mean()
        d_ge[i] = baseline_ge - global_efficiency(A, [i])

    # Node strength is the control that decides whether any of this is new.
    # Gu et al. (2015) report average controllability correlating strongly with
    # weighted degree. If the resection ranking is just strength, then nothing
    # here beats summing a row of the connectivity matrix.
    strength = A.sum(axis=1)

    return {"d_ac": d_ac, "d_mc": d_mc, "d_ge": d_ge, "pagerank": pr,
            "strength": strength}


def rank_desc(values):
    order = np.argsort(-np.asarray(values, dtype=float))
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(1, len(order) + 1)
    return ranks


def partial_spearman(x, y, control):
    """Spearman correlation of x and y after regressing out `control`.

    Answers the question that matters: does controllability tell you anything
    about resection damage that node strength does not already tell you?
    """
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rc = np.argsort(np.argsort(control)).astype(float)
    design = np.column_stack([np.ones_like(rc), rc])
    res_x = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    res_y = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(np.corrcoef(res_x, res_y)[0, 1])


def report(res, label, top_frac=0.1):
    d_ac, d_ge, pr = res["d_ac"], res["d_ge"], res["pagerank"]
    strength = res["strength"]
    n = len(d_ac)
    cutoff = max(1, int(round(n * top_frac)))

    rho_ac, p_ac = spearmanr(d_ac, d_ge)
    rho_mc, _ = spearmanr(res["d_mc"], d_ge)
    rho_pr, _ = spearmanr(pr, d_ge)

    rho_ac_str, _ = spearmanr(d_ac, strength)
    rho_ge_str, _ = spearmanr(d_ge, strength)
    partial = partial_spearman(d_ac, d_ge, strength)

    r_ac, r_ge = rank_desc(d_ac), rank_desc(d_ge)
    only_ac = np.flatnonzero((r_ac <= cutoff) & (r_ge > cutoff))
    only_ge = np.flatnonzero((r_ge <= cutoff) & (r_ac > cutoff))

    worst_ge = int(np.argmax(d_ge))
    worst_ac = int(np.argmax(d_ac))
    pr_order = np.argsort(-pr)

    print(f"=== {label} ===")
    print(f"  nodes {n}, edge density "
          f"{np.count_nonzero(np.triu(res['A'], 1)) / (n * (n - 1) / 2):.3f}")
    print(f"  Spearman(d_ave_control, d_global_eff) = {rho_ac:+.3f} (p={p_ac:.2g})")
    print(f"  Spearman(d_mod_control, d_global_eff) = {rho_mc:+.3f}")
    print(f"  Spearman(pagerank,      d_global_eff) = {rho_pr:+.3f}")
    print(f"  -- redundancy check against node strength --")
    print(f"  Spearman(d_ave_control, strength)     = {rho_ac_str:+.3f}"
          f"   {'REDUNDANT' if abs(rho_ac_str) > 0.95 else ''}")
    print(f"  Spearman(d_global_eff,  strength)     = {rho_ge_str:+.3f}")
    print(f"  partial Spearman(d_ac, d_ge | strength) = {partial:+.3f}")
    print(f"  worst parcel by GE: {worst_ge} | by controllability: {worst_ac} "
          f"| same: {worst_ge == worst_ac}")
    print(f"  highest-PageRank parcel is worst by GE: "
          f"top-1 {pr_order[0] == worst_ge}, top-3 {worst_ge in pr_order[:3]}")
    print(f"  top {int(top_frac * 100)}% disagreement ({cutoff} per tier): "
          f"{len(only_ac)} flagged by controllability only, "
          f"{len(only_ge)} by GE only")
    print(f"    controllability-only parcels: {only_ac.tolist()[:12]}"
          f"{' ...' if len(only_ac) > 12 else ''}")
    print()


def main():
    for label, thr in [("dense (no threshold)", None),
                       ("top 15% of edges", 0.15),
                       ("top 5% of edges", 0.05)]:
        A = load_rosen_halgren(CSV, threshold=thr)
        res = sweep_single_parcel(A)
        res["A"] = A
        report(res, label)


if __name__ == "__main__":
    main()
