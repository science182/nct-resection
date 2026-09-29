"""Kept so `python3 audit.py ...` still works. The tool lives in nct_resection.audit,
and once installed it is also available as the `nct-audit` command."""

from nct_resection.audit import *  # noqa: F401,F403
from nct_resection.audit import main

if __name__ == "__main__":
    main()
