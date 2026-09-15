#!/usr/bin/env python3
"""Execute both paper notebooks and export HTML copies with embedded figures."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    for builder in ["build_paper_notebook.py", "build_appendix_notebook.py"]:
        subprocess.run([sys.executable, str(ROOT / "scripts" / builder), "--execute"], cwd=ROOT, check=True)
    subprocess.run([
        sys.executable, "-m", "nbconvert", "--to", "html", "--embed-images",
        "--output-dir", str(ROOT / "results/notebooks"),
        str(ROOT / "for_rtm3_paper.ipynb"), str(ROOT / "for_appendix_RTM3.ipynb"),
    ], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
