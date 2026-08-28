"""Filesystem roots, resolved from the environment.

No machine-specific path is hard-coded anywhere in this repository. Set these once:

    export VAAD_ROOT=/path/to/your/data-root     # required for anything touching raw data
    export VAAD_TOOLS=$VAAD_ROOT/tools           # third-party models/binaries (optional)

`scripts/00_download_dependencies.sh` creates the expected layout under $VAAD_ROOT.
"""
import os

VAAD_ROOT = os.environ.get("VAAD_ROOT", os.path.expanduser("~/vaad-data"))
VAAD_TOOLS = os.environ.get("VAAD_TOOLS", os.path.join(VAAD_ROOT, "tools"))

AACDB_STRUCTURES = os.path.join(VAAD_ROOT, "datasets", "AACDB", "complex_structure")
RESULTS_DIR = os.path.join(VAAD_ROOT, "results")
WORK_DIR = os.path.join(VAAD_ROOT, "work")


def require(path: str, what: str) -> str:
    """Fail loudly and early when an external dependency is not staged."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{what} not found at {path}. Set VAAD_ROOT/VAAD_TOOLS and run "
            f"scripts/00_download_dependencies.sh (see data/README.md)."
        )
    return path
