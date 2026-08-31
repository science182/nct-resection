"""Two confound checks the spin test does not cover.

1. Parcel geometry. The spin test rules out spatial autocorrelation but not
   size. Larger parcels contain more vertices and absorb more streamlines, so
   deleting one perturbs the network differently for reasons that have nothing
   to do with position in the network. If the top-scoring parcels are
   systematically large or small, geometry explains the result.

2. Residualization choice. Degree was removed with a quadratic fit on
   within-subject strength ranks. That was a judgement call. If the language
   enrichment appears under that choice and vanishes under a linear fit or a
   plain rank correction, the finding is an artifact of the correction rather
   than of the data.
"""

import numpy as np
from scipy.stats import spearmanr

from annot import read_annot, read_gifti_surface
from run_energy import load_network_assignment
from spin_test import parcel_centroids, spin_permutations

K = 36


def parcel_geometry():
    """Vertex count and surface area per parcel, in connectome row order."""
    verts = np.zeros(360)
    areas = np.zeros(360)
    for offset, hemi, surf in ((0, "lh", "L"), (180, "rh", "R")):
        labels, _ = read_annot(f"data/{hemi}.HCP-MMP1.annot")
        coords, faces = read_gifti_surface(f"data/fsaverage_{surf}_white.surf.gii")

        tri = coords[faces]
        tri_area = 0.5 * np.linalg.norm(
            np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        # Assign each triangle to the parcel holding the majority of its vertices.
        tri_lab = labels[faces]
        owner = np.where(tri_lab[:, 0] == tri_lab[:, 1], tri_lab[:, 0], tri_lab[:, 2])

        for p in range(1, 181):
            verts[offset + p - 1] = np.sum(labels == p)
            areas[offset + p - 1] = tri_area[owner == p].sum()
    return verts, areas


def residualize(damage, strength, method):
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


def language_enrichment(score, asg, perms):
    members = set(np.flatnonzero((asg == 5) | (asg == 15)).tolist())
    obs = len(set(np.argsort(-score)[:K].tolist()) & members)
    spun = np.array([len(set(np.argsort(-score[p])[:K].tolist()) & members)
                     for p in perms])
    return obs, (spun >= obs).mean()


def main():
    d = np.load("data/subject_deletions.npz")
    d_ac, st = d["d_ac"], d["strength"]
    asg, _ = load_network_assignment()

    print("=== 1. parcel geometry ===")
    verts, areas = parcel_geometry()
    res = residualize(d_ac, st, "quadratic").mean(0)
    top = np.argsort(-res)[:K]
    rest = np.setdiff1d(np.arange(360), top)

    for name, g in (("vertex count", verts), ("surface area mm2", areas)):
        rho = spearmanr(res, g)[0]
        print(f"  {name:<18} corr with residual {rho:+.3f} | "
              f"top{K} median {np.median(g[top]):8.0f} vs rest {np.median(g[rest]):8.0f} "
              f"(ratio {np.median(g[top]) / np.median(g[rest]):.2f})")
    print(f"  strength corr with surface area: {spearmanr(st.mean(0), areas)[0]:+.3f}")

    verdict = ("geometry does not explain the ranking"
               if abs(spearmanr(res, areas)[0]) < 0.25
               else "GEOMETRY CONFOUND, size tracks the score")
    print(f"  -> {verdict}\n")

    print("=== 2. does the language result depend on how degree was removed? ===")
    print("  building spin permutations...")
    perms = spin_permutations(parcel_centroids("lh", "L"),
                              parcel_centroids("rh", "R"),
                              n_spins=1000)
    print(f"  {'method':<14} {'Language in top36':>18} {'p spin':>9} "
          f"{'resid vs strength':>19}")
    for method in ("linear", "quadratic", "cubic", "rank-only"):
        r = residualize(d_ac, st, method)
        prof = r.mean(0)
        obs, p = language_enrichment(prof, asg, perms)
        leak = np.mean([spearmanr(r[s], st[s])[0] for s in range(r.shape[0])])
        print(f"  {method:<14} {obs:>13} of {K} {p:>9.4f} {leak:>19.3f}")
    print("\n  expected by chance: 2.3 of 36")


if __name__ == "__main__":
    main()
