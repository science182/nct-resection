"""Does the language finding replicate on raw streamline counts?

Three questions, in order of how badly a negative answer would hurt:

1. Does degree-corrected controllability still concentrate on language cortex
   when the connectome is weighted by raw streamline counts rather than Fpt?
2. Does it hold in two independent halves of the subjects? This is the
   replication the Fpt result could not provide at n=24.
3. Does the stronger size confound in raw counts explain it? Raw counts are not
   normalized per seed, so parcel size feeds directly into strength.
"""

import numpy as np
from scipy.stats import spearmanr

from run_confounds import K, language_enrichment, parcel_geometry, residualize
from run_energy import load_network_assignment
from spin_test import parcel_centroids, spin_permutations

FPT = "data/subject_deletions.npz"
STREAM = "data/streamline_deletions.npz"


def report(tag, d_ac, st, asg, perms, methods=("quadratic",)):
    for method in methods:
        r = residualize(d_ac, st, method)
        prof = r.mean(0)
        obs, p = language_enrichment(prof, asg, perms)
        leak = np.mean([spearmanr(r[s], st[s])[0] for s in range(r.shape[0])])
        print(f"  {tag:<28} {method:<10} {obs:>2} of {K}  p_spin {p:.4f}  "
              f"resid~strength {leak:+.3f}")
    return residualize(d_ac, st, "quadratic")


def main():
    asg, _ = load_network_assignment()
    print("building spin permutations...")
    perms = spin_permutations(parcel_centroids("lh", "L"),
                              parcel_centroids("rh", "R"), n_spins=1000)

    fpt = np.load(FPT)
    stream = np.load(STREAM)
    n_stream = stream["d_ac"].shape[0]
    print(f"Fpt: {fpt['d_ac'].shape[0]} subjects | "
          f"streamline: {n_stream} subjects\n")

    print("=== 1. does it replicate on streamline counts? (expect 2.3 of 36 by chance)")
    res_fpt = report("Fpt (original)", fpt["d_ac"], fpt["strength"], asg, perms)
    res_str = report("streamline counts", stream["d_ac"], stream["strength"],
                     asg, perms, methods=("linear", "quadratic", "rank-only"))

    print("\n=== 2. do the two weightings agree on which parcels are risky?")
    rho = spearmanr(res_fpt.mean(0), res_str.mean(0))[0]
    print(f"  Spearman(Fpt residual map, streamline residual map) = {rho:+.3f}")

    print("\n=== 3. independent split of the streamline subjects")
    half = n_stream // 2
    a = report(f"subjects 1-{half}", stream["d_ac"][:half],
               stream["strength"][:half], asg, perms)
    b = report(f"subjects {half + 1}-{n_stream}", stream["d_ac"][half:],
               stream["strength"][half:], asg, perms)
    print(f"  map agreement between halves: "
          f"{spearmanr(a.mean(0), b.mean(0))[0]:+.3f}")

    print("\n=== 4. size confound, which raw counts should worsen")
    verts, areas = parcel_geometry()
    for name, res, st in (("Fpt", res_fpt, fpt["strength"]),
                          ("streamline", res_str, stream["strength"])):
        prof = res.mean(0)
        top = np.argsort(-prof)[:K]
        rest = np.setdiff1d(np.arange(360), top)
        print(f"  {name:<11} strength~area {spearmanr(st.mean(0), areas)[0]:+.3f} | "
              f"residual~area {spearmanr(prof, areas)[0]:+.3f} | "
              f"top{K} size ratio {np.median(areas[top]) / np.median(areas[rest]):.2f}")


if __name__ == "__main__":
    main()
