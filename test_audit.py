"""Tests for audit.py. Run with `python3 test_audit.py`. No pytest required.

Both bugs this file guards against shipped in the first version of the tool, and
both were found only because it was run on a second dataset:

  The .mat loader tested `dtype == object` before checking for structured
  dtypes. A dtype like [('SC','O'),('FC','O')] is not equal to object, so every
  MATLAB struct was silently skipped and the loader reported no matrix at all.

  The recommendation suggested log10(1 + w) unconditionally. That only
  compresses weights above 1. On a connectome whose weights are all below 1 it
  is nearly linear, and it moved a tail ratio of 193 to 184 while claiming to
  fix the problem.

The second is the more dangerous kind: it produced confident, plausible,
useless advice rather than an error.
"""

import os
import sys
import tempfile

import numpy as np
import scipy.io as sio

from audit import (
    _from_mat,
    _parcel_count,
    _square_numeric,
    clean,
    damage_sweep,
    load_any,
    power_transform,
    recommend_alpha,
    tail_ratio,
)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def heavy_tailed(n=60, seed=0):
    """Connectome with a wide, skewed weight distribution."""
    rng = np.random.default_rng(seed)
    W = rng.lognormal(0.0, 3.0, size=(n, n))
    W = np.triu(W, 1)
    return clean(W + W.T)


def test_structured_mat_is_found():
    """The exact layout that defeated the first loader."""
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "nested.mat")
        inner = np.empty((2, 1), dtype=object)
        inner[0, 0] = heavy_tailed(30)
        inner[1, 0] = heavy_tailed(40)
        outer = np.empty((1, 1), dtype=[("SC", "O"), ("FC", "O")])
        outer[0, 0]["SC"] = inner
        outer[0, 0]["FC"] = np.zeros((5, 5))
        sio.savemat(path, {"connMatrices": outer, "labels": np.array(["a", "b"])})

        found = _from_mat(path, max_nodes=400)
        check("matrix recovered from a nested MATLAB struct",
              found is not None and _parcel_count(found.shape) in (30, 40),
              f"got shape {None if found is None else found.shape}")


def test_prefers_scale_under_cap():
    """Multi-scale files must not silently select the largest parcellation."""
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "scales.mat")
        arr = np.empty((3, 1), dtype=object)
        for i, n in enumerate((60, 200, 900)):
            arr[i, 0] = heavy_tailed(n)
        sio.savemat(path, {"SC": arr})

        found = _from_mat(path, max_nodes=400)
        check("largest scale at or under the cap is chosen",
              _parcel_count(found.shape) == 200,
              f"chose {_parcel_count(found.shape)}")

        found_small = _from_mat(path, max_nodes=100)
        check("cap is respected when lowered",
              _parcel_count(found_small.shape) == 60,
              f"chose {_parcel_count(found_small.shape)}")


def test_label_arrays_are_not_mistaken_for_connectomes():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "withlabels.mat")
        sio.savemat(path, {"parcelIDs": np.array([["L_V1"], ["L_MST"]], dtype=object),
                           "conn": heavy_tailed(25)})
        found = _from_mat(path, max_nodes=400)
        check("string labels ignored, connectome found",
              found is not None and _parcel_count(found.shape) == 25)


def test_square_numeric_rejects_non_numeric():
    check("object dtype rejected",
          not _square_numeric(np.empty((4, 4), dtype=object)))
    check("string dtype rejected",
          not _square_numeric(np.array([["a", "b"], ["c", "d"]])))
    check("float square accepted", _square_numeric(np.zeros((6, 6))))
    check("non-square rejected", not _square_numeric(np.zeros((3, 7))))


def test_load_any_orients_stacks():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "stack.npy")
        np.save(path, np.stack([heavy_tailed(20) for _ in range(4)]))
        A = load_any(path)
        check("subject-major stack preserved", A.shape == (4, 20, 20),
              f"got {A.shape}")

        path2 = os.path.join(d, "single.npy")
        np.save(path2, heavy_tailed(20))
        B = load_any(path2)
        check("single matrix gains a subject axis", B.shape == (1, 20, 20),
              f"got {B.shape}")


def test_log10_advice_would_have_been_caught():
    """The scale-dependence bug, stated as a test.

    log10(1 + w) compresses weights far above 1 and does nothing to weights far
    below it. A recommendation engine that ignores this gives useless advice on
    normalized connectomes.
    """
    # The failing regime is when EVERY weight is below 1, as in a connectome
    # normalized by streamline length and region size, where the Lausanne
    # release tops out at about 0.08. A fixture that straddles 1 does not test
    # this: log10(1 + w) still compresses whatever part lies above 1.
    raw = heavy_tailed(40)
    big = raw / raw.max() * 1e6
    small = raw / raw.max() * 0.08
    assert small.max() < 1.0

    before_big, after_big = tail_ratio(big), tail_ratio(np.log10(1.0 + big))
    before_small, after_small = tail_ratio(small), tail_ratio(np.log10(1.0 + small))

    check("log10(1+w) compresses weights above 1",
          after_big < before_big / 10,
          f"{before_big:.0f} -> {after_big:.0f}")
    check("log10(1+w) barely moves weights below 1 (the bug)",
          after_small > before_small / 2,
          f"{before_small:.0f} -> {after_small:.0f}")


def test_power_transform_is_scale_free():
    """The fix: a power transform must work at either scale."""
    for label, scale in (("large", 1e6), ("small", 1e-3)):
        M = heavy_tailed(40) * scale
        alpha = recommend_alpha([M], target_tail=30.0)
        ok = alpha is not None
        achieved = tail_ratio(power_transform(M, alpha)) if ok else float("inf")
        check(f"power transform reaches the target on {label} weights",
              ok and achieved < 30.0,
              f"alpha={alpha}, tail {tail_ratio(M):.0f} -> {achieved:.1f}")


def test_power_transform_preserves_edges_and_order():
    M = heavy_tailed(30)
    P = power_transform(M, 0.5)
    check("no edges added or removed", np.array_equal(M > 0, P > 0))
    off = ~np.eye(30, dtype=bool)
    a, b = M[off], P[off]
    keep = a > 0
    check("weight ordering preserved",
          np.array_equal(np.argsort(a[keep]), np.argsort(b[keep])))


def test_tail_ratio_direction():
    M = heavy_tailed(40)
    check("compression lowers the tail ratio",
          tail_ratio(power_transform(M, 0.25)) < tail_ratio(M))
    uniform = clean(np.ones((20, 20)))
    check("a flat connectome has tail ratio 1",
          abs(tail_ratio(uniform) - 1.0) < 1e-9,
          f"got {tail_ratio(uniform)}")


def test_damage_sweep_subsamples():
    M = heavy_tailed(50)
    full, st_full = damage_sweep(M, "ac")
    sub, st_sub = damage_sweep(M, "ac", max_parcels=10)
    check("full sweep covers every parcel", full.shape == (50,), f"{full.shape}")
    check("subsampled sweep honours the cap", sub.shape == (10,), f"{sub.shape}")
    check("strength vector matches the sweep length", st_sub.shape == sub.shape)
    check("damage is positive for a connected graph", bool(np.all(full > 0)))


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


# --- resection_report checks -------------------------------------------------

def test_ties_are_not_wins():
    """Identical candidates must grade IDENTICAL, not ROBUST.

    Growing from different seeds can converge on the same corridor. Counting
    only strict wins made every tie a win for whichever candidate came second,
    so two identical resections were reported as a reliable difference.
    """
    from resection_report import grade_pairs
    scores = {("raw", "strength"): np.array([1.0, 1.0]),
              ("log", "strength"): np.array([1.0, 1.0])}
    row = grade_pairs(scores, [[1], [1]])[0]
    check("identical candidates grade IDENTICAL", row["grade"] == "IDENTICAL",
          f"got {row['grade']} at {row['share']:.0%}")


def test_grades_track_agreement():
    from resection_report import grade_pairs
    unanimous = {f"c{k}": np.array([1.0, 2.0]) for k in range(10)}
    check("unanimous comparison is ROBUST",
          grade_pairs(unanimous, [[1], [2]])[0]["grade"] == "ROBUST")

    split = {f"c{k}": np.array([1.0, 2.0]) if k < 5 else np.array([2.0, 1.0])
             for k in range(10)}
    check("evenly split comparison is a COIN FLIP",
          grade_pairs(split, [[1], [2]])[0]["grade"] == "COIN FLIP")

    mostly = {f"c{k}": np.array([1.0, 2.0]) if k < 9 else np.array([2.0, 1.0])
              for k in range(10)}
    row = grade_pairs(mostly, [[1], [2]])[0]
    check("nine of ten is LIKELY", row["grade"] == "LIKELY",
          f"got {row['grade']} at {row['share']:.0%}")


def test_damage_is_monotone_in_extent():
    """Removing more tissue cannot be scored as less damaging."""
    from resection_report import damage
    M = heavy_tailed(40)
    small = damage(M, [3], "strength")
    large = damage(M, [3, 4, 5], "strength")
    check("strength damage grows with resection size", large > small)
    e_small = damage(M, [3], "efficiency")
    e_large = damage(M, [3, 4, 5], "efficiency")
    check("efficiency damage grows with resection size", e_large > e_small,
          f"{e_small:.4g} vs {e_large:.4g}")

if __name__ == "__main__":
    sys.exit(main())
