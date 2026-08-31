"""Minimum control energy for state transitions, and how a resection changes it.

Average and modal controllability say how easy a node is to drive in general.
Control energy asks a sharper question: to move the brain from one specific
state to another specific state, how much input does it take, and how much more
does it take after tissue is removed.

Discrete-time minimum energy for x0 -> xT in T steps:

    W_T = sum_{k=0}^{T-1} A^k B B' (A')^k
    E   = (xT - A^T x0)' W_T^+ (xT - A^T x0)
"""

import numpy as np

from controllability import normalize_adjacency, spectral_scale


def controllability_gramian(A, B, horizon):
    """Finite-horizon controllability Gramian. A must already be normalized."""
    n = A.shape[0]
    W = np.zeros((n, n))
    Ak = np.eye(n)
    BBt = B @ B.T
    for _ in range(horizon):
        W += Ak @ BBt @ Ak.T
        Ak = Ak @ A
    return W


def minimum_control_energy(A, x0, xT, B=None, horizon=10, normalize=True,
                           ref_scale=None):
    """Energy to drive the network from x0 to xT in `horizon` steps.

    B defaults to the identity, meaning every node can receive input. Pass a
    diagonal 0/1 matrix to restrict control to a subset of nodes.

    The Gramian is often ill-conditioned, so the linear system is solved in the
    least-squares sense, which is equivalent to using the pseudo-inverse.
    """
    A = (normalize_adjacency(A, ref_scale=ref_scale) if normalize
         else np.asarray(A, dtype=float))
    n = A.shape[0]
    B = np.eye(n) if B is None else np.asarray(B, dtype=float)
    x0 = np.asarray(x0, dtype=float).ravel()
    xT = np.asarray(xT, dtype=float).ravel()
    if x0.size != n or xT.size != n:
        raise ValueError(f"state vectors must have length {n}")

    W = controllability_gramian(A, B, horizon)
    delta = xT - np.linalg.matrix_power(A, horizon) @ x0
    z, *_ = np.linalg.lstsq(W, delta, rcond=None)
    return float(delta @ z)


def apply_resection(A, resected, mode="disconnect"):
    """Return the post-resection connectivity matrix.

    Two ways to simulate removing tissue, and the choice matters:

    "disconnect" zeroes the resected node's edges but keeps it in the matrix.
    The network stays the same size, so state vectors and energies before and
    after are directly comparable. This is the right choice for control energy.

    "delete" removes the rows and columns entirely. Closer to how graph metrics
    are usually reported, but it changes the dimension, which makes energies
    between the two networks not strictly comparable.
    """
    A = np.asarray(A, dtype=float)
    resected = np.asarray(sorted(set(resected)), dtype=int)

    if mode == "disconnect":
        out = A.copy()
        out[resected, :] = 0.0
        out[:, resected] = 0.0
        return out
    if mode == "delete":
        retained = np.setdiff1d(np.arange(A.shape[0]), resected)
        return A[np.ix_(retained, retained)]
    raise ValueError(f"unknown mode {mode!r}")


def delta_control_energy(A, resected, x0, xT, horizon=10, B=None):
    """Extra energy required to reach the same target state after resection.

    Uses "disconnect" so both networks have the same dimension and the same
    state vectors, and normalizes both by the intact network's scale for the
    same reason delta_controllability does. A positive `extra` means the
    resection made the transition more expensive.
    """
    A = np.asarray(A, dtype=float)
    ref = spectral_scale(A)

    intact = minimum_control_energy(A, x0, xT, B=B, horizon=horizon,
                                    ref_scale=ref)
    lesioned_A = apply_resection(A, resected, mode="disconnect")
    lesioned = minimum_control_energy(lesioned_A, x0, xT, B=B, horizon=horizon,
                                      ref_scale=ref)

    return {
        "intact": intact,
        "lesioned": lesioned,
        "extra": lesioned - intact,
        "ratio": lesioned / intact if intact else float("nan"),
    }


def energies_for_targets(A, x0, targets, horizon=10, B=None, normalize=True,
                         ref_scale=None):
    """Minimum energy for many target states at once.

    The Gramian depends only on the network and the horizon, not on the target,
    so it is built once and reused. That turns a per-resection sweep over ten
    target states from ten matrix-power loops into one.
    """
    A = (normalize_adjacency(A, ref_scale=ref_scale) if normalize
         else np.asarray(A, dtype=float))
    n = A.shape[0]
    B = np.eye(n) if B is None else np.asarray(B, dtype=float)

    W = controllability_gramian(A, B, horizon)
    decayed = np.linalg.matrix_power(A, horizon) @ np.asarray(x0, dtype=float)

    out = []
    for xT in targets:
        delta = np.asarray(xT, dtype=float).ravel() - decayed
        z, *_ = np.linalg.lstsq(W, delta, rcond=None)
        out.append(float(delta @ z))
    return np.asarray(out)


def network_state(assignments, active_networks, n=None):
    """Build a target state vector from network membership.

    `assignments` is a length-N array naming each node's network. Nodes in
    `active_networks` get 1, everything else 0. This is a placeholder for a real
    clinically defined state, and which states are worth scoring is exactly the
    open question worth asking a clinician rather than deciding here.
    """
    assignments = np.asarray(assignments)
    n = assignments.size if n is None else n
    active = set(active_networks)
    return np.array([1.0 if a in active else 0.0 for a in assignments[:n]])
