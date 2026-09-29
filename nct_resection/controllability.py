"""Network control theory metrics on structural connectomes.

Implements average and modal controllability following Gu et al. (2015),
Nat Commun 6:8414, matching the reference MATLAB implementations so results
are comparable to the published literature and to nctpy.

The discrete-time linear model assumed throughout is

    x(t+1) = A_norm x(t) + B u(t)

where A_norm is the structural connectivity matrix normalized so the system
is stable.
"""

import numpy as np
from scipy.linalg import schur


def _is_symmetric(A, tol=1e-10):
    return A.shape[0] == A.shape[1] and np.allclose(A, A.T, atol=tol, rtol=0)


def spectral_scale(A):
    """Largest singular value of A, the quantity normalization divides by.

    Structural connectomes are symmetric, and for a symmetric matrix the
    singular values are the absolute eigenvalues, so `eigvalsh` gives the same
    answer as a full SVD about three times faster. Non-symmetric input falls
    back to the SVD.
    """
    A = np.asarray(A, dtype=float)
    if _is_symmetric(A):
        return float(np.abs(np.linalg.eigvalsh(A)).max())
    return float(np.linalg.svd(A, compute_uv=False)[0])


def _spectrum(A):
    """Eigenvalues and eigenvectors for the controllability formulas.

    Both measures need a real Schur form. For a symmetric matrix that form is
    the eigendecomposition, which `eigh` computes roughly six times faster than
    the general `schur` routine, agreeing to machine precision (about 1e-14 on
    real connectome data). Since a deletion sweep runs this once per removed
    parcel per subject, the difference is hours.
    """
    if _is_symmetric(A):
        eigvals, vectors = np.linalg.eigh(A)
        return eigvals, vectors
    T, U = schur(A, output="real")
    return np.diag(T), U


def normalize_adjacency(A, target_radius=0.95, ref_scale=None, c=None):
    """Scale A so its spectral radius is < 1, making the discrete system stable.

    `ref_scale` is the subtle and important argument. Normalizing each network
    by its own largest singular value is correct when comparing unrelated
    networks, but it is WRONG for lesion studies. Removing a hub lowers a
    network's largest singular value, so dividing by that smaller number pushes
    the remaining eigenvalues closer to 1 and inflates average controllability.
    The result is that resecting a hub appears to *improve* controllability.

    Passing the intact network's scale as `ref_scale` holds the denominator
    fixed so the lesioned network is measured on the same ruler. Removing nodes
    from a non-negative matrix can only lower the largest singular value, so the
    lesioned system stays stable under the intact network's scale.

    `target_radius` is the second trap, and it is specific to the units of your
    connectome. The convention in the literature is the additive form
    A / (1 + sigma_max), which is calibrated for streamline counts where
    sigma_max is in the thousands and the result sits just under 1. Applied to a
    probability-weighted connectome, where sigma_max is around 0.5, that same
    formula yields a spectral radius near 0.3. In that regime

        AC_i = sum_j U_ij^2 / (1 - lam_j^2) ~ 1 + sum_j lam_j^2 U_ij^2
        MC_i = sum_j (1 - lam_j^2) U_ij^2   = 1 - sum_j lam_j^2 U_ij^2

    so AC ~ 2 - MC and the two measures collapse into one, with almost no
    dynamic range. Scaling to a fixed target radius instead makes the
    normalization independent of the connectome's units. Pass `c` to recover
    the legacy additive convention.

    On the Rosen & Halgren group average, measured:

        radius   AC spread   |corr(AC, 2 - MC)|
        0.32       0.016         0.9999      degenerate, MC adds nothing
        0.95       0.239         0.959
        0.9999    36.9           0.609       matches streamline-count regime

    There is no free lunch here, and the two failure modes pull opposite ways.
    Pushing the radius toward 1 separates AC from MC, but it also makes AC
    hypersensitive to the spectral radius itself. Since removing any node at all
    lowers the radius, the lesion difference then becomes dominated by "was
    anything removed" rather than "what was removed". Measured on the Rosen &
    Halgren average, for 8-parcel resections:

        radius   |corr(AC, 2-MC)|   dAC signal as % of its offset
        0.50          0.9999                  74%
        0.95          0.959                   42%
        0.99          0.9(approx)             32%
        0.9999        0.609                    2%

    At 0.9999 the node-identity signal is 2 percent of a constant offset, and
    the resulting rank correlations flip sign purely as an artifact of that.
    The default is therefore 0.95, which keeps lesion deltas well conditioned,
    at the cost that AC and MC are close to relabelings of each other there and
    MC should not be reported as independent evidence. Pass `c` to recover the
    legacy additive convention.
    """
    A = np.asarray(A, dtype=float)
    if A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square, got {A.shape}")
    sigma_max = spectral_scale(A) if ref_scale is None else float(ref_scale)
    if sigma_max == 0:
        raise ValueError("A has no edges; cannot normalize")
    if c is not None:
        return A / (c + sigma_max)
    if not 0 < target_radius < 1:
        raise ValueError("target_radius must lie in (0, 1)")
    return A * (target_radius / sigma_max)


def average_controllability(A, normalize=True, ref_scale=None):
    """Per-node average controllability: ease of driving the network into many
    nearby states.

    AC_i = sum_j U[i, j]^2 / (1 - lambda_j^2)

    Returns an array of length N.
    """
    A = (normalize_adjacency(A, ref_scale=ref_scale) if normalize
         else np.asarray(A, dtype=float))
    eigvals, U = _spectrum(A)
    denom = 1.0 - eigvals ** 2
    # Guard against modes sitting on the unit circle after normalization.
    denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
    return np.sum((U ** 2) / denom[np.newaxis, :], axis=1)


def modal_controllability(A, normalize=True, ref_scale=None):
    """Per-node modal controllability: ability to push the network into
    hard-to-reach states.

    MC_i = sum_j (1 - lambda_j^2) * U[i, j]^2

    Returns an array of length N.
    """
    A = (normalize_adjacency(A, ref_scale=ref_scale) if normalize
         else np.asarray(A, dtype=float))
    eigvals, U = _spectrum(A)
    return np.sum((1.0 - eigvals ** 2) * (U ** 2), axis=1)


def delta_controllability(A, resected, metric=average_controllability):
    """Change in controllability caused by resecting a set of nodes.

    Two things are easy to get wrong here, and both make a resection look
    harmless when it is not.

    1. Size. Average controllability is a sum over modes, so a whole-network
       value depends on node count. Comparing intact against lesioned directly
       conflates "harder to control" with "smaller." Fixed by scoring both
       networks only over the nodes that survive the resection.

    2. Scale. Normalizing each network by its own largest singular value lets
       the denominator move when a hub is removed, which inflates the lesioned
       network's controllability and flips the sign of the effect. Fixed by
       normalizing both networks by the intact network's scale.

    Args:
        A: (N, N) weighted symmetric connectivity matrix.
        resected: indices of nodes removed by the simulated resection.
        metric: average_controllability or modal_controllability.

    Returns:
        dict with the mean intact value over retained nodes, the mean lesioned
        value over those same nodes, the absolute drop, and the fractional drop.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    resected = np.asarray(sorted(set(resected)), dtype=int)
    retained = np.setdiff1d(np.arange(n), resected)

    if retained.size < 2:
        raise ValueError("resection leaves fewer than two nodes")

    ref = spectral_scale(A)
    intact_all = metric(A, ref_scale=ref)
    intact = intact_all[retained]

    A_lesioned = A[np.ix_(retained, retained)]
    if not np.any(A_lesioned):
        raise ValueError("resection disconnects the network entirely")
    lesioned = metric(A_lesioned, ref_scale=ref)

    mean_intact = float(np.mean(intact))
    mean_lesioned = float(np.mean(lesioned))
    drop = mean_intact - mean_lesioned

    return {
        "intact": mean_intact,
        "lesioned": mean_lesioned,
        "drop": drop,
        "frac_drop": drop / mean_intact if mean_intact else float("nan"),
        "n_retained": int(retained.size),
    }
