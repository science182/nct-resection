"""Does control energy escape the degree baseline where controllability did not?

Average controllability turned out to track node strength at +0.77 to +0.91
regardless of resection size, which makes it a degree proxy. Control energy is
the one remaining quantity in this project that is not a spectral summary of
the connectivity matrix: it depends on a specific target state, so two parcels
with identical strength can differ if one sits on the path to that state and
the other does not.

Target states are bilateral activation of each of the ten networks shipped with
the Rosen & Halgren ordering file. Resections are simulated by disconnecting
rather than deleting, so the state vector keeps its length and energies before
and after are comparable.
"""

import h5py
import numpy as np
from scipy.stats import spearmanr

from controllability import spectral_scale
from data import load_rosen_halgren
from energy import apply_resection, energies_for_targets
from lesion import global_efficiency
from run_real import CSV, partial_spearman

ORDER_MAT = "data/parcelOrder_and_networkAssignment.mat"
HORIZON = 10


def load_network_assignment(path=ORDER_MAT):
    """Per-parcel network index in the connectome's own row order.

    The file stores the network index in reordered space alongside the
    permutation, so it has to be inverted rather than used directly.
    """
    with h5py.File(path, "r") as f:
        labels = [
            "".join(chr(c) for c in np.array(f[r]).ravel())
            for r in np.array(f["networkLabels"]).ravel()
        ]
        reordered_idx = np.array(f["postReorderingNetworkIndex"]).ravel().astype(int)
        order = np.array(f["networkOrder"]).ravel().astype(int)

    assignment = np.zeros(360, dtype=int)
    assignment[order - 1] = reordered_idx
    return assignment, [s.strip() for s in labels]


def bilateral_targets(assignment, labels):
    """One target state per network, both hemispheres active."""
    names, states = [], []
    for i in range(10):
        left, right = i + 1, i + 11
        mask = (assignment == left) | (assignment == right)
        names.append(labels[i][2:].strip())
        states.append(mask.astype(float))
    return names, states


def main():
    A = load_rosen_halgren(CSV)
    assignment, labels = load_network_assignment()
    names, targets = bilateral_targets(assignment, labels)
    n = A.shape[0]
    ref = spectral_scale(A)
    x0 = np.zeros(n)

    print("parcels per bilateral network:")
    print("  " + ", ".join(f"{nm} {int(t.sum())}" for nm, t in zip(names, targets)))

    intact = energies_for_targets(A, x0, targets, horizon=HORIZON, ref_scale=ref)
    print(f"\nintact energies: min {intact.min():.4g}, max {intact.max():.4g}")
    if not np.all(np.isfinite(intact)):
        raise SystemExit("non-finite intact energies; aborting")

    strength = A.sum(axis=1)
    baseline_ge = global_efficiency(A)

    extra = np.zeros((n, len(targets)))
    d_ge = np.zeros(n)
    for i in range(n):
        lesioned = apply_resection(A, [i], mode="disconnect")
        extra[i] = energies_for_targets(lesioned, x0, targets, horizon=HORIZON,
                                        ref_scale=ref) - intact
        d_ge[i] = baseline_ge - global_efficiency(A, [i])

    mean_extra = extra.mean(axis=1)

    print("\nper-network: does extra control energy just track node strength?")
    print(f"  {'network':<18} {'rho vs strength':>16} {'rho vs dGE':>12} "
          f"{'excl. in-network':>18}")
    for j, nm in enumerate(names):
        outside = assignment != (j + 1)
        outside &= assignment != (j + 11)
        rho_s = spearmanr(extra[:, j], strength)[0]
        rho_g = spearmanr(extra[:, j], d_ge)[0]
        rho_out = spearmanr(extra[outside, j], strength[outside])[0]
        print(f"  {nm:<18} {rho_s:>+16.3f} {rho_g:>+12.3f} {rho_out:>+18.3f}")

    print("\naveraged over all ten networks:")
    print(f"  Spearman(extra energy, strength)      = "
          f"{spearmanr(mean_extra, strength)[0]:+.3f}")
    print(f"  Spearman(extra energy, d_global_eff)  = "
          f"{spearmanr(mean_extra, d_ge)[0]:+.3f}")
    print(f"  partial (extra energy, d_ge | strength) = "
          f"{partial_spearman(mean_extra, d_ge, strength):+.3f}")

    print("\nis the ranking target-specific, or the same parcels every time?")
    rhos = [spearmanr(extra[:, a], extra[:, b])[0]
            for a in range(len(names)) for b in range(a + 1, len(names))]
    print(f"  median cross-network rank agreement: {np.median(rhos):+.3f}")
    print(f"  range: {min(rhos):+.3f} to {max(rhos):+.3f}")


if __name__ == "__main__":
    main()
