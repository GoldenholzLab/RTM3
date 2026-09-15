#!/usr/bin/env python3
"""Build the methodological appendix using shared reporting functions."""
import argparse
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]


def build_notebook():
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    cells = [md((ROOT / "appendix_methods.md").read_text()),
        code("from pathlib import Path\nimport matplotlib.pyplot as plt\nfrom rtm3_reporting import *\nresults = load_results(Path('results'))\nshow_design(results)"),
        md("## Primary results\n\nEvery condition is reported. RTM means the observed trajectory indicator. Positive median percentage change, MPC, indicates apparent improvement without treatment."),
        code("show_outcomes(results['primary'])"),
        md("## Selection diagnostics\n\nFrequencies are seizures per 30-day month among participants meeting primary criteria."),
        code("show_selection(results)"),
        md("## Fixed and reselected cohorts\n\nThe fixed cohort is selected once under perfect detection. MPC applies only to positive-baseline participants; zero baselines are reported separately."),
        code("show_outcomes(results['cohort'], 'cohort_strategy')\ncomparison_figure(results['cohort'], 'cohort_strategy', 'Fixed versus reselected cohorts')\nplt.show()"),
        md("## Additive expected-alarm thresholds\n\nOnly eligibility uses raw counts. Outcomes retain monthly correction."),
        code("threshold_comparison(results)"),
        md("## The 25-day rule\n\nDifferences are sensitivity analysis minus main analysis, in percentage points."),
        code("show_outcomes(results['interval'], 'eligibility_variant')\ninterval_comparison(results)\nplt.show()"),
        md("## Correction and interval alternatives\n\nThe raw-outcome comparison includes expected alarms in its reference. Oracle removal uses simulated alarm identities."),
        code("show_outcomes(results['correction'], 'correction_strategy')\ncomparison_figure(results['correction'], 'correction_strategy', 'False-alarm correction and interval alternatives')\nplt.show()"),
        md("## Interpretation"), code("summarize_results(results)"),
        md("## Reproducibility checks"), code("validate_artifacts(results)")]
    return nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    notebook = build_notebook()
    if args.execute:
        NotebookClient(notebook, timeout=1200, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}}).execute()
    output = ROOT / "for_appendix_RTM3.ipynb"
    nbf.write(notebook, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
