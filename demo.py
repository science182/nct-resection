"""End-to-end run on a synthetic connectome.

The connectome here is FAKE. It is a modular small-world graph with roughly
brain-like degree structure, and it exists only so the pipeline runs today and
you can see the shape of the output. Every number this prints is meaningless
biologically. Replacing `synthetic_connectome` with real HCP-MMP1 connectivity
matrices is the entire remaining scientific step.
"""

import numpy as np

from compare import summarize, sweep_resections


def synthetic_connectome(n=180, n_modules=6, seed=0):
    """Modular weighted symmetric graph, standing in for one hemisphere."""
    rng = np.random.default_rng(seed)
    module = np.repeat(np.arange(n_modules), int(np.ceil(n / n_modules)))[:n]

    same = module[:, None] == module[None, :]
    p = np.where(same, 0.35, 0.03)
    edges = rng.random((n, n)) < p
    edges = np.triu(edges, k=1)
    edges = edges | edges.T

    weights = rng.lognormal(mean=0.0, sigma=0.6, size=(n, n))
    weights = np.triu(weights, k=1)
    weights = weights + weights.T

    A = np.where(edges, weights, 0.0)

    # A few hubs, so PageRank has something to find.
    hubs = rng.choice(n, size=8, replace=False)
    for h in hubs:
        targets = rng.choice(n, size=n // 6, replace=False)
        A[h, targets] = A[targets, h] = rng.lognormal(1.2, 0.4, size=targets.size)
    np.fill_diagonal(A, 0.0)
    return A


def main():
    A = synthetic_connectome()
    print(f"synthetic connectome: {A.shape[0]} nodes, "
          f"{int((A > 0).sum() / 2)} edges\n")

    for size in (1, 4, 8):
        print(f"--- contiguous resections of {size} parcel(s) ---")
        rows = sweep_resections(A, size=size)
        print(summarize(rows))
        print()


if __name__ == "__main__":
    main()
