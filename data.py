"""Loading real connectomes.

The best public match for this project is Rosen & Halgren (2021), eNeuro
8(1):ENEURO.0416-20.2020, "A Whole-Cortex Probabilistic Diffusion Tractography
Connectome". It is HCP-MMP1.0, 360 cortical parcels, 1065 HCP subjects, and it
publishes both the group average and every individual matrix under CC-BY 4.0.

    https://doi.org/10.5281/zenodo.4060485

Files worth knowing about in that record:

    averageConnectivity_Fpt.csv          1.0 MB   group average, plain CSV
    individualConnectivity_10^Fpt.mat    1.0 GB   1065 subjects, MATLAB v7.3
    parcelOrder_and_networkAssignment.mat        parcel order + network labels

Start with the CSV. The individual file is what eventually makes the
interindividual variability question answerable, and MATLAB v7.3 is HDF5
underneath, so h5py reads it without MATLAB.

One thing to verify before trusting any result: "Fpt" is a fractional
probability of tractography, and the individual file is named 10^Fpt, which
implies the stored values may be log-scaled. Confirm whether the CSV holds Fpt
or 10^Fpt before using the weights, because controllability is sensitive to the
weight distribution and a silent log/linear mixup changes every ranking.
"""

import numpy as np

ZENODO_DOI = "https://doi.org/10.5281/zenodo.4060485"
N_PARCELS_HCP_MMP1 = 360


def load_matrix(path, delimiter=",", skip_header=0, skip_cols=0):
    """Load a dense connectivity matrix from a text/CSV file."""
    M = np.genfromtxt(path, delimiter=delimiter, skip_header=skip_header)
    if skip_cols:
        M = M[:, skip_cols:]
    return np.asarray(M, dtype=float)


def load_rosen_halgren(path, threshold=None):
    """Load averageConnectivity_Fpt.csv and return usable linear weights.

    Verified against the real file: 360x360, symmetric, NaN on the diagonal,
    and the stored values are log10 fractional probabilities running from about
    -6.2 to -0.8. They MUST be exponentiated. Controllability requires
    non-negative weights, and feeding the raw log values in would either error
    out or, worse, silently invert the weight ordering.

    After exponentiating, the matrix is fully dense: every parcel pair carries
    some tractography probability, spanning six orders of magnitude. `threshold`
    keeps only the strongest fraction of edges (e.g. 0.1 for the top 10 percent),
    which is standard for probabilistic connectomes and worth varying, since
    every downstream ranking depends on it.
    """
    M = np.genfromtxt(path, delimiter=",")
    if M.shape[0] != M.shape[1]:
        raise ValueError(f"expected square matrix, got {M.shape}")
    if np.nanmax(M) > 0:
        raise ValueError(
            "expected log10 values (all negative); this file may already be "
            "linear, in which case skip the exponentiation"
        )

    A = np.power(10.0, M)
    A[~np.isfinite(A)] = 0.0
    np.fill_diagonal(A, 0.0)
    A = (A + A.T) / 2.0

    if threshold is not None:
        A = threshold_edges(A, threshold)
    return A


def threshold_edges(A, keep_fraction):
    """Keep the strongest `keep_fraction` of edges, zero the rest."""
    A = np.asarray(A, dtype=float).copy()
    iu = np.triu_indices_from(A, k=1)
    w = A[iu]
    nonzero = w[w > 0]
    if nonzero.size == 0:
        return A
    cutoff = np.quantile(nonzero, 1.0 - keep_fraction)
    mask = A < cutoff
    A[mask] = 0.0
    np.fill_diagonal(A, 0.0)
    return (A + A.T) / 2.0


def prepare(A, symmetrize=True, zero_diagonal=True, drop_isolated=True):
    """Put a raw matrix into the shape the rest of this code assumes.

    Returns the cleaned matrix and the indices of nodes kept, so parcel labels
    can be realigned if any node is dropped.
    """
    A = np.asarray(A, dtype=float).copy()
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"expected a square matrix, got {A.shape}")

    A[~np.isfinite(A)] = 0.0
    if A.min() < 0:
        raise ValueError("negative edge weights; controllability assumes "
                         "non-negative connectivity")
    if symmetrize:
        A = (A + A.T) / 2.0
    if zero_diagonal:
        np.fill_diagonal(A, 0.0)

    kept = np.arange(A.shape[0])
    if drop_isolated:
        connected = A.sum(axis=1) > 0
        if not connected.all():
            kept = np.flatnonzero(connected)
            A = A[np.ix_(kept, kept)]
    return A, kept


def describe(A):
    """Quick sanity summary. Run this before trusting anything downstream."""
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    offdiag = A[~np.eye(n, dtype=bool)]
    nonzero = offdiag[offdiag > 0]
    density = nonzero.size / offdiag.size if offdiag.size else 0.0
    return {
        "nodes": n,
        "density": density,
        "symmetric": bool(np.allclose(A, A.T)),
        "weight_min": float(nonzero.min()) if nonzero.size else 0.0,
        "weight_max": float(nonzero.max()) if nonzero.size else 0.0,
        "weight_median": float(np.median(nonzero)) if nonzero.size else 0.0,
    }


def load_individual_connectomes(path, key=None):
    """Read the 1065-subject MATLAB v7.3 file. Requires h5py.

    Returns an array shaped (subjects, parcels, parcels). The file is 1 GB, so
    load it once and cache what you need rather than re-reading per analysis.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "h5py is required to read MATLAB v7.3 files: pip install h5py"
        ) from exc

    with h5py.File(path, "r") as f:
        keys = [k for k in f.keys() if not k.startswith("#")]
        if key is None:
            if len(keys) != 1:
                raise ValueError(f"specify key; file contains {keys}")
            key = keys[0]
        arr = np.array(f[key])

    # MATLAB writes column-major, so HDF5 gives back transposed axes.
    if arr.ndim == 3 and arr.shape[0] == arr.shape[1]:
        arr = np.transpose(arr, (2, 0, 1))
    return arr
