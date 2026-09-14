"""Does the instability extend beyond control measures to standard graph metrics?

Everything so far concerns average controllability and global efficiency in a
deletion design. If weighting choice only destabilized those, the finding would
be narrow. This tests a panel of the nodal measures that structural connectomics
actually reports, asking the question papers actually ask: which regions does
this analysis flag as important?

Nine measures, all standard:

    strength                total connectivity, the trivial baseline
    degree                  binary connection count
    betweenness             fraction of shortest paths passing through
    closeness               inverse mean distance to all others
    eigenvector             recursive centrality
    pagerank                damped recursive centrality
    clustering              weighted local triangle density
    communicability         sum over weighted walks of all lengths
    avg controllability     ease of driving the network to nearby states
    modal controllability   ease of driving it to distant states

Each is computed once per weighting convention, which is far cheaper than a
deletion sweep, and is also closer to what is reported: hub identification and
centrality rankings appear more often in this literature than deletion damage.

For every measure we ask how much its ranking changes when only the edge
weighting changes, using the same statement types as run_fliprate.py.

    python3 run_measures.py
"""

import itertools
import os

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.stats import spearmanr

from controllability import average_controllability, modal_controllability
from run_weightings import WEIGHTINGS, reweight

TOP_N = 10


def _distance(A):
    with np.errstate(divide="ignore"):
        D = np.where(A > 0, 1.0 / A, np.inf)
    np.fill_diagonal(D, 0.0)
    return D


def nodal_measures(A):
    """Every measure as a length-N vector. Higher means more central."""
    n = A.shape[0]
    G = nx.from_numpy_array(A)
    out = {}

    out["strength"] = A.sum(axis=1)
    out["degree"] = (A > 0).sum(axis=1).astype(float)

    # Shortest-path measures share one distance solve.
    sp = shortest_path(_distance(A), method="FW", directed=False)
    finite = np.isfinite(sp)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(finite & (sp > 0), 1.0 / sp, 0.0)
    out["closeness"] = inv.sum(axis=1) / max(n - 1, 1)

    bt = nx.betweenness_centrality(G, weight="distance_proxy", normalized=True) \
        if False else nx.betweenness_centrality(
            nx.from_numpy_array(_distance_for_nx(A)), weight="weight",
            normalized=True)
    out["betweenness"] = np.array([bt[i] for i in range(n)])

    try:
        ev = nx.eigenvector_centrality_numpy(G, weight="weight")
        out["eigenvector"] = np.array([ev[i] for i in range(n)])
    except Exception:
        out["eigenvector"] = np.zeros(n)

    pr = nx.pagerank(G, weight="weight")
    out["pagerank"] = np.array([pr[i] for i in range(n)])

    cl = nx.clustering(G, weight="weight")
    out["clustering"] = np.array([cl[i] for i in range(n)])

    # Communicability: sum of weighted walks, on the spectrally normalized
    # matrix so the exponential does not overflow on streamline counts.
    scale = np.abs(np.linalg.eigvalsh(A)).max()
    out["communicability"] = expm(A / scale).sum(axis=1) if scale > 0 \
        else np.zeros(n)

    out["avg_control"] = average_controllability(A)
    out["mod_control"] = modal_controllability(A)
    return out


def _distance_for_nx(A):
    """Finite distance matrix for betweenness; unreachable pairs get no edge."""
    D = _distance(A)
    D[~np.isfinite(D)] = 0.0
    return D


def jaccard(a, b):
    return len(set(a.tolist()) & set(b.tolist())) / \
           len(set(a.tolist()) | set(b.tolist()))


def compare(measures_by_weighting, names, conventions=WEIGHTINGS):
    """Per measure, stability of its ranking across weighting conventions."""
    rows = []
    for m in names:
        pairs = list(itertools.combinations(conventions, 2))
        rhos, jacs, hubs = [], [], []
        for a, b in pairs:
            x = measures_by_weighting[a][m]
            y = measures_by_weighting[b][m]
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            rhos.append(spearmanr(x, y)[0])
            jacs.append(jaccard(np.argsort(-x)[:TOP_N], np.argsort(-y)[:TOP_N]))
            hubs.append(int(np.argmax(x) == np.argmax(y)))
        if not rhos:
            continue
        rows.append({
            "measure": m,
            "rho": float(np.mean(rhos)),
            "top10": float(np.mean(jacs)),
            "hub_changes": 1.0 - float(np.mean(hubs)),
        })
    return rows


WEIGHTED_ONLY = ("raw", "log", "rownorm")


def run(label, loader, n_subjects=4, conventions=WEIGHTINGS):
    """Average ranking stability over several subjects.

    `conventions` defaults to all four. Passing WEIGHTED_ONLY drops binarizing,
    which is the most violent change and could be said to do the work on its
    own. The finding survives either way.
    """
    acc = {}
    for i in range(n_subjects):
        A = loader(i)
        measures = {w: nodal_measures(reweight(A, w)) for w in conventions}
        rows = compare(measures, list(measures[conventions[0]]),
                       conventions=conventions)
        for r in rows:
            acc.setdefault(r["measure"], []).append(r)

    summary = [{
        "measure": m,
        "rho": float(np.mean([r["rho"] for r in rs])),
        "top10": float(np.mean([r["top10"] for r in rs])),
        "hub_changes": float(np.mean([r["hub_changes"] for r in rs])),
    } for m, rs in acc.items()]

    print(f"\n=== {label} ===")
    print(f"  {n_subjects} subjects, conventions: {', '.join(conventions)}")
    print(f"  {'measure':<18} {'mean rho':>9} {'top-10 overlap':>15} "
          f"{'top hub changes':>16}")
    for r in sorted(summary, key=lambda r: -r["rho"]):
        print(f"  {r['measure']:<18} {r['rho']:>+9.3f} {r['top10']:>15.0%} "
              f"{r['hub_changes']:>16.0%}")
    low = sum(1 for r in summary if r["top10"] < 0.7)
    print(f"  measures whose top-10 overlap is below 70%: {low} of {len(summary)}")
    # Binary degree must be invariant under monotone reweighting, since raw,
    # log and rownorm change weights without changing which edges exist. It is
    # NOT expected to hold once binarizing is included, because that convention
    # thresholds edges away, so degree legitimately changes. Checking it there
    # would be a false alarm.
    if "binary" not in conventions and any(r["measure"] == "degree"
                                           for r in summary):
        d = next(r for r in summary if r["measure"] == "degree")
        ok = d["rho"] > 0.999 and d["top10"] > 0.999
        print(f"  internal check, binary degree invariant under monotone "
              f"reweighting: {'yes' if ok else 'NO, pipeline problem'}")
    return summary


def main():
    from audit import clean, load_any

    hcp_all = None
    lau_all = None

    def hcp(i):
        nonlocal hcp_all
        if hcp_all is None:
            hcp_all = load_any("data/streamlineCount.mat")
        return clean(hcp_all[i])

    def lausanne(i):
        nonlocal lau_all
        if lau_all is None:
            lau_all = load_any("data/lausanne70.mat", max_nodes=400)
        return clean(lau_all[i])

    print("\nDOES WEIGHTING DESTABILIZE STANDARD GRAPH METRICS TOO?")
    print("  nodal centrality rankings, which is what hub papers report")

    if os.path.exists("data/streamlineCount.mat"):
        run("HCP streamline counts, all four conventions", hcp)
        run("HCP streamline counts, weighted conventions only", hcp,
            conventions=WEIGHTED_ONLY)
    else:
        print("\n  skipped HCP: data/streamlineCount.mat not present")

    if os.path.exists("data/lausanne70.mat"):
        run("Lausanne 219, all four conventions", lausanne)
    else:
        print("\n  skipped Lausanne: data/lausanne70.mat not present")

    print("\n  A high mean rho with a low top-10 overlap is the common case:")
    print("  the bulk ordering survives while the specific regions a paper")
    print("  would name do not. Binary degree is included as an internal check,")
    print("  since it must be invariant to any monotone reweighting.")


if __name__ == "__main__":
    main()
