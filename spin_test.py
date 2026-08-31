"""Spin test: is the language enrichment real, or just spatial autocorrelation?

The permutation null used earlier shuffled network labels freely. That is
anti-conservative for brain maps. The controllability residual is spatially
smooth and language parcels are spatially clustered, so a smooth map and
clustered labels will produce apparent enrichment under free shuffling even
when nothing is there.

The spin test (Alexander-Bloch et al. 2018; parcel-level variant after
Vasa et al. 2018) fixes this by rotating the map on the cortical sphere. Random
rotations preserve the map's spatial autocorrelation and the parcellation's
geometry, and only destroy the alignment between the map and the anatomy. If
enrichment survives, the alignment is doing the work rather than the smoothness.

Rotated parcels are matched to originals with the Hungarian algorithm so the
result is a genuine one-to-one permutation, rather than the cheaper
nearest-neighbour matching that lets one parcel absorb several others.
"""

import numpy as np
from scipy.optimize import linear_sum_assignment

from annot import read_annot, read_gifti_surface


def parcel_centroids(hemi, surf):
    """Centroid of each parcel on the sphere, projected back to the surface."""
    labels, _ = read_annot(f"data/{hemi}.HCP-MMP1.annot")
    coords, _ = read_gifti_surface(f"data/fsaverage_{surf}_sphere.surf.gii")

    cents = np.zeros((180, 3))
    for p in range(1, 181):
        pts = coords[labels == p]
        c = pts.mean(0)
        cents[p - 1] = c / np.linalg.norm(c) * 100.0
    return cents


def random_rotation(rng):
    """Uniform random 3D rotation via QR, with the reflection removed."""
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def spin_permutations(cents_l, cents_r, n_spins=1000, seed=0):
    """One-to-one parcel permutations induced by random sphere rotations.

    The same rotation is applied to both hemispheres, mirrored in x for the
    right, which keeps left and right treated consistently.
    """
    rng = np.random.default_rng(seed)
    flip = np.diag([-1.0, 1.0, 1.0])
    perms = np.zeros((n_spins, 360), dtype=int)

    for s in range(n_spins):
        R = random_rotation(rng)
        for offset, cents, rot in ((0, cents_l, R),
                                   (180, cents_r, flip @ R @ flip)):
            rotated = cents @ rot.T
            cost = np.linalg.norm(rotated[:, None, :] - cents[None, :, :], axis=2)
            _, col = linear_sum_assignment(cost)
            perms[s, offset:offset + 180] = col + offset
    return perms


def main():
    import h5py

    from run_energy import load_network_assignment
    from run_position import residualize_within_subject

    d = np.load("data/subject_deletions.npz")
    res_ac = residualize_within_subject(d["d_ac"], d["strength"]).mean(0)
    res_ge = residualize_within_subject(d["d_ge"], d["strength"]).mean(0)
    asg, labels = load_network_assignment()

    print("building parcel centroids on the sphere...")
    cents_l = parcel_centroids("lh", "L")
    cents_r = parcel_centroids("rh", "R")

    n_spins = 1000
    print(f"generating {n_spins} spin permutations (Hungarian matching)...")
    perms = spin_permutations(cents_l, cents_r, n_spins=n_spins)

    uniq = np.mean([len(np.unique(p)) for p in perms])
    print(f"  mean distinct parcels per permutation: {uniq:.1f} / 360 "
          f"(360 means genuinely one-to-one)\n")

    rng = np.random.default_rng(1)
    K = 36
    for score_name, score in (("controllability residual", res_ac),
                              ("global efficiency residual", res_ge)):
        print(f"=== {score_name}, top {K} parcels ===")
        top = set(np.argsort(-score)[:K].tolist())
        print(f"  {'network':<20} {'obs':>4} {'exp':>5} {'p naive':>9} {'p spin':>8}")
        for net in range(1, 11):
            members = np.flatnonzero((asg == net) | (asg == net + 10))
            obs = len(top & set(members.tolist()))
            if obs < 2:
                continue

            # Naive null: free label shuffling, ignores spatial structure.
            naive = np.array([
                len(set(rng.choice(360, K, replace=False)) & set(members.tolist()))
                for _ in range(2000)
            ])
            # Spin null: rotate the map, preserving spatial autocorrelation.
            spun = np.array([
                len(set(np.argsort(-score[p])[:K].tolist()) & set(members.tolist()))
                for p in perms
            ])
            p_naive = (naive >= obs).mean()
            p_spin = (spun >= obs).mean()
            flag = "" if p_spin < 0.05 else "  <- does not survive"
            print(f"  {labels[net - 1].strip()[2:].strip():<20} {obs:>4} "
                  f"{K * len(members) / 360:>5.1f} {p_naive:>9.4f} "
                  f"{p_spin:>8.4f}{flag}")
        print()


if __name__ == "__main__":
    main()
