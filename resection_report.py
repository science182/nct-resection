"""Kept so `python3 resection_report.py ...` still works. The tool lives in
nct_resection.report, and once installed it is also available as the
`nct-resection-report` command."""

from nct_resection.report import *  # noqa: F401,F403
from nct_resection.report import main

if __name__ == "__main__":
    main()
