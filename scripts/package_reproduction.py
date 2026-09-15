#!/usr/bin/env python3
"""Copy the canonical analysis into a standalone submission reproduction folder."""
import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = ["README.md", "LICENSE", "requirements.txt", "requirements-reproduction.txt",
              "reproduce_rtm3.py", "realSim.py", "rtm3_reporting.py", "appendix_methods.md",
              "for_rtm3_paper.ipynb", "for_appendix_RTM3.ipynb", "old_for_rtm3_paper.ipynb"]
SCRIPT_FILES = ["build_paper_notebook.py", "build_appendix_notebook.py", "reproduce_all.py"]


def package(destination):
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ROOT_FILES:
        shutil.copy2(ROOT / filename, destination / filename)
    for directory, filenames in [("scripts", SCRIPT_FILES), ("tests", ["test_reproduce_rtm3.py"])]:
        (destination / directory).mkdir(exist_ok=True)
        for filename in filenames:
            shutil.copy2(ROOT / directory / filename, destination / directory / filename)
    output = destination / "results"
    output.mkdir(exist_ok=True)
    for path in (ROOT / "results").iterdir():
        if path.is_file() and (path.name.startswith("rtm3_") or path.name.startswith("figure") or path.name == "appendix_selection_diagnostics.png"):
            shutil.copy2(path, output / path.name)
    if (ROOT / "results/notebooks").exists():
        shutil.copytree(ROOT / "results/notebooks", output / "notebooks", dirs_exist_ok=True)
    print(f"Packaged canonical reproduction in {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    package(parser.parse_args().destination.resolve())
