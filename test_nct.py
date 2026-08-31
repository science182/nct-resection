"""Tests. Run with `python3 test_nct.py`. No pytest required.

The first three tests exist because of a real bug. An earlier version of this
code normalized each lesioned network by its own largest singular value, which
made resecting a hub look protective and produced a confident, clean, entirely
inverted result. Nothing caught it. These tests catch it.
"""

import sys

import numpy as np

from controllability import (
    average_controllability,
    delta_controllability,
    modal_controllability,
    normalize_adjacency,
    spectral_scale,
)
from demo import synthetic_connectome
from energy import (
    apply_resection,
    controllability_gramian,
    delta_control_energy,
    minimum_control_energy,
)
from lesion import grow_resection, pagerank_scores

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def hub_and_leaf(A):
    pr = pagerank_scores(A)
    return int(np.argmax(pr)), int(np.argmin(pr))


def test_normalization_stabilizes():
    A = synthetic_connectome(n=60, seed=1)
    An = normalize_adjacency(A)
    radius = max(abs(np.linalg.eigvals(An)))
    check("normalized spectral radius < 1", radius < 1.0, f"got {radius:.4f}")


def test_hub_costs_more_than_leaf():
    """The headline sanity check. Removing a hub must cost more controllability
    than removing a peripheral node. The old normalization bug inverted this."""
    A = synthetic_connectome(n=120, seed=2)
    hub, leaf = hub_and_leaf(A)
    hub_drop = delta_controllability(A, [hub])["drop"]
    leaf_drop = delta_controllability(A, [leaf])["drop"]
    check(
        "hub removal costs more average controllability than leaf removal",
        hub_drop > leaf_drop,
        f"hub {hub_drop:+.5f} vs leaf {leaf_drop:+.5f}",
    )
    check(
        "removing a hub reduces controllability at all",
        hub_drop > 0,
        f"drop was {hub_drop:+.5f}",
    )


def test_per_network_normalization_reproduces_the_bug():
    """Documents the failure mode so it cannot quietly come back.

    Normalizing each network by its own scale should visibly understate, or
    invert, the cost of removing a hub compared to the shared-scale answer.
    """
    A = synthetic_connectome(n=120, seed=2)
    hub, _ = hub_and_leaf(A)
    retained = np.setdiff1d(np.arange(A.shape[0]), [hub])
    sub = A[np.ix_(retained, retained)]

    shared = spectral_scale(A)
    correct = (average_controllability(A, ref_scale=shared)[retained].mean()
               - average_controllability(sub, ref_scale=shared).mean())
    buggy = (average_controllability(A)[retained].mean()
             - average_controllability(sub).mean())

    check(
        "per-network normalization understates hub cost (the old bug)",
        buggy < correct,
        f"buggy {buggy:+.5f} vs correct {correct:+.5f}",
    )


def test_ac_and_mc_are_not_degenerate():
    """Average and modal controllability must carry distinct information.

    If the spectral radius after normalization is too far below 1, then
    AC ~ 1 + s and MC = 1 - s for the same s, so AC ~ 2 - MC and modal
    controllability becomes a relabeling of average controllability. This
    happened for real on a probability-weighted connectome under the legacy
    additive normalization.
    """
    A = synthetic_connectome(n=120, seed=11)
    ac = average_controllability(A)
    mc = modal_controllability(A)
    corr = abs(np.corrcoef(ac, 2 - mc)[0, 1])
    check("AC is not a relabeling of MC", corr < 0.99, f"|corr(AC, 2-MC)| = {corr:.5f}")

    legacy_ac = average_controllability(A, ref_scale=None)
    check("default normalization gives usable dynamic range",
          (ac.max() - ac.min()) > 1e-2,
          f"spread {ac.max() - ac.min():.2e}")


def test_delta_signal_is_not_swamped_by_offset():
    """The lesion difference must reflect which nodes were removed.

    Near marginal stability, removing any node at all drops the spectral radius
    and collapses average controllability, so every resection produces nearly
    the same delta. The node-identity signal then rides as a few percent on a
    large constant, and rank correlations computed on it are noise that can flip
    sign. Guarded because this produced a clean, strong, entirely spurious
    result on real data.
    """
    A = synthetic_connectome(n=100, seed=12)
    drops = np.array([delta_controllability(A, [i])["drop"] for i in range(0, 100, 5)])
    signal = (drops.max() - drops.min()) / abs(drops.mean())
    check("delta carries node-specific signal, not just an offset",
          signal > 0.10, f"spread is {signal:.1%} of the offset")


def test_delta_compares_same_nodes():
    A = synthetic_connectome(n=80, seed=3)
    resected = grow_resection(A, seed=5, size=6)
    result = delta_controllability(A, resected)
    check(
        "delta scored over retained nodes only",
        result["n_retained"] == 80 - len(resected),
        f"got {result['n_retained']}, expected {80 - len(resected)}",
    )


def test_modal_controllability_range():
    A = synthetic_connectome(n=60, seed=4)
    mc = modal_controllability(A)
    check("modal controllability in [0, 1]",
          bool(np.all(mc >= 0) and np.all(mc <= 1)),
          f"range [{mc.min():.4f}, {mc.max():.4f}]")


def test_gramian_positive_semidefinite():
    A = normalize_adjacency(synthetic_connectome(n=50, seed=5))
    W = controllability_gramian(A, np.eye(50), horizon=8)
    eigs = np.linalg.eigvalsh((W + W.T) / 2)
    check("Gramian is positive semidefinite", eigs.min() > -1e-9,
          f"min eigenvalue {eigs.min():.3e}")


def test_energy_non_negative():
    A = synthetic_connectome(n=50, seed=6)
    rng = np.random.default_rng(0)
    x0 = np.zeros(50)
    xT = rng.random(50)
    e = minimum_control_energy(A, x0, xT, horizon=8)
    check("minimum control energy is non-negative", e >= 0, f"got {e:.4f}")


def test_disconnect_preserves_shape():
    A = synthetic_connectome(n=40, seed=7)
    out = apply_resection(A, [3, 9], mode="disconnect")
    deleted = apply_resection(A, [3, 9], mode="delete")
    check("disconnect keeps matrix dimension", out.shape == A.shape,
          f"got {out.shape}")
    check("disconnect zeroes the resected edges",
          bool(np.all(out[3] == 0) and np.all(out[:, 9] == 0)))
    check("delete shrinks matrix dimension", deleted.shape == (38, 38),
          f"got {deleted.shape}")


def test_hub_resection_costs_more_energy():
    A = synthetic_connectome(n=90, seed=8)
    hub, leaf = hub_and_leaf(A)
    rng = np.random.default_rng(1)
    x0 = np.zeros(90)
    xT = rng.random(90)

    hub_extra = delta_control_energy(A, [hub], x0, xT, horizon=8)["extra"]
    leaf_extra = delta_control_energy(A, [leaf], x0, xT, horizon=8)["extra"]
    check(
        "resecting a hub raises control energy more than resecting a leaf",
        hub_extra > leaf_extra,
        f"hub {hub_extra:+.5f} vs leaf {leaf_extra:+.5f}",
    )


def test_resection_is_contiguous():
    A = synthetic_connectome(n=100, seed=9)
    resected = grow_resection(A, seed=12, size=5)
    adjacency = A > 0
    # every node after the first must touch an earlier one
    ok = all(adjacency[node, resected[:i]].any()
             for i, node in enumerate(resected) if i > 0)
    check("grown resection is connected", ok, f"got {resected}")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    print(f"running {len(tests)} test groups\n")
    for t in tests:
        print(f"{t.__name__}:")
        t()
        print()

    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {', '.join(FAILURES)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
