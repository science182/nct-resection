"""Edge-weighting conventions and degree correction.

These are the two operations every result in this project depends on, so they
live in the package rather than in the analysis scripts that first used them.
`run_weightings.py` and `run_confounds.py` re-export them unchanged, so older
imports keep working.
"""

import numpy as np

WEIGHTINGS = ("raw", "log", "rownorm", "binary")


def reweight(M, kind):
    """M arrives symmetrized with a zero diagonal."""
    if kind == "raw":
        return M
    if kind == "log":
        return np.log10(1.0 + M)
    if kind == "rownorm":
        rs = M.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        R = M / rs
        return (R + R.T) / 2.0
    if kind == "binary":
        iu = np.triu_indices_from(M, k=1)
        w = M[iu]
        nz = w[w > 0]
        cut = np.quantile(nz, 0.70)
        B = (M >= cut).astype(float)
        np.fill_diagonal(B, 0.0)
        return B
    raise ValueError(kind)


def residualize(damage, strength, method):
    """Remove the within-subject strength-rank component from damage ranks.

    `damage` and `strength` are (subjects, parcels). `method` is one of
    "linear", "quadratic", "cubic" or "rank-only".
    """
    out = np.zeros_like(damage)
    for s in range(damage.shape[0]):
        x = np.argsort(np.argsort(strength[s])).astype(float)
        y = np.argsort(np.argsort(damage[s])).astype(float)
        if method == "rank-only":
            out[s] = y - x
            continue
        cols = [np.ones_like(x), x]
        if method == "quadratic":
            cols.append(x ** 2)
        elif method == "cubic":
            cols.extend([x ** 2, x ** 3])
        design = np.column_stack(cols)
        out[s] = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return out
