"""Fetch every input this project needs. Nothing here is redistributed.

    python3 download_data.py           # core files, about 13 MB
    python3 download_data.py --full    # adds the per-subject connectomes, 1.3 GB

Core is enough to run the tests, the demo, and the group-average analyses.
The full set is only needed for the per-subject and streamline-weighting work.

Sources:
  Rosen & Halgren (2021), eNeuro 8(1), CC-BY 4.0
      https://doi.org/10.5281/zenodo.4060485
  fsaverage surfaces via TemplateFlow
      https://github.com/templateflow/tpl-fsaverage
"""

import argparse
import os
import sys
import urllib.request

ZENODO = "https://zenodo.org/records/4060485/files/{}?download=1"
TEMPLATEFLOW = ("https://templateflow.s3.amazonaws.com/tpl-fsaverage/"
                "tpl-fsaverage_hemi-{}_den-164k_{}.surf.gii")

# (destination filename, url, approx MB, needed for)
CORE = [
    ("averageConnectivity_Fpt.csv", ZENODO.format("averageConnectivity_Fpt.csv"),
     1.0, "group-average connectome, log10 Fpt"),
    ("parcelOrder_and_networkAssignment.mat",
     ZENODO.format("parcelOrder_and_networkAssignment.mat"),
     0.2, "parcel order and the ten networks"),
    ("lh.HCP-MMP1.annot", ZENODO.format("lh.HCP-MMP1.annot"),
     1.3, "left parcellation, for labels and adjacency"),
    ("rh.HCP-MMP1.annot", ZENODO.format("rh.HCP-MMP1.annot"),
     1.3, "right parcellation"),
    ("fsaverage_L_white.surf.gii", TEMPLATEFLOW.format("L", "white"),
     4.3, "left surface faces, for parcel adjacency"),
    ("fsaverage_R_white.surf.gii", TEMPLATEFLOW.format("R", "white"),
     4.3, "right surface faces"),
    ("fsaverage_L_sphere.surf.gii", TEMPLATEFLOW.format("L", "sphere"),
     3.3, "left sphere, for the spin test"),
    ("fsaverage_R_sphere.surf.gii", TEMPLATEFLOW.format("R", "sphere"),
     3.3, "right sphere"),
]

FULL = [
    ("individualConnectivity.mat",
     ZENODO.format("individualConnectivity_10%5EFpt.mat"),
     1000.0, "1065 per-subject connectomes, linear 10^Fpt"),
    ("streamlineCount.mat",
     ZENODO.format("individualConnectivity_rawStreamlineCount.mat"),
     308.0, "1065 per-subject raw streamline counts"),
]


def fetch(dest_dir, name, url, mb, purpose):
    path = os.path.join(dest_dir, name)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print(f"  have    {name}")
        return True
    print(f"  getting {name}  (~{mb:.0f} MB, {purpose})", flush=True)
    tmp = path + ".part"
    try:
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, path)
        return True
    except Exception as exc:  # noqa: BLE001 - report and continue
        if os.path.exists(tmp):
            os.remove(tmp)
        print(f"  FAILED  {name}: {exc}")
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="also fetch the 1.3 GB per-subject connectomes")
    ap.add_argument("--dir", default="data")
    args = ap.parse_args()

    os.makedirs(args.dir, exist_ok=True)
    targets = CORE + (FULL if args.full else [])
    total = sum(t[2] for t in targets)
    print(f"fetching {len(targets)} files into {args.dir}/ (~{total:.0f} MB)\n")

    failed = [name for name, url, mb, why in targets
              if not fetch(args.dir, name, url, mb, why)]

    print()
    if failed:
        print(f"{len(failed)} file(s) failed: {', '.join(failed)}")
        return 1
    print("all files present")
    if not args.full:
        print("run with --full for the per-subject analyses "
              "(run_subjects, run_subject_deletions, run_streamline, run_weightings)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
