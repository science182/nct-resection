"""Node-deletion analysis of structural connectomes, with the checks that decide
whether a deletion-based map can be trusted.

Library modules:
    controllability  average and modal controllability (Gu et al. 2015)
    lesion           contiguous resection growth, weighted global efficiency, PageRank
    energy           minimum control energy and resection application
    compare          rank resections by every metric and extract disagreements
    weighting        edge-weighting conventions and degree correction
    annot            FreeSurfer .annot / GIFTI readers and HCP-MMP1 parcel adjacency
    audit            per-dataset audit (also the `nct-audit` command)
    report           convention-robust comparison of candidate resections
                     (also the `nct-resection-report` command)

The analysis scripts that produced the results in FINDINGS.md live at the top
level of the repository, not in this package.
"""

__version__ = "0.1.0"
