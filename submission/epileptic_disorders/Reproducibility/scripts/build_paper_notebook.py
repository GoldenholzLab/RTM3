#!/usr/bin/env python3
"""Build the function-driven reproduction notebook for both manuscript figures."""
import argparse
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]


def build_notebook():
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    cells = [md("# RTM3 main-paper reproduction\n\nThis notebook reruns the 100,000-participant simulation with fixed seeds, generates both manuscript images through shared functions, and embeds the figures and source results. Figure 1 is calculated from outcomes; Figure 2 is conceptual. The Appendix contains equations and sensitivity analyses. The archived `old_for_rtm3_paper.ipynb` is not a source of current results."),
        code("from pathlib import Path\nimport matplotlib.pyplot as plt\nfrom IPython.display import Image, display\nfrom reproduce_rtm3 import run_analysis, make_main_figure, DEFAULT_SEED\nfrom rtm3_reporting import load_results, show_outcomes, show_main_cohort_comparison, make_conceptual_figure, validate_artifacts\nOUTPUT = Path('results')"),
        code("paths = run_analysis(n_patients=100_000, seed=DEFAULT_SEED, output_dir=OUTPUT, progress_every=25_000)\nresults = load_results(OUTPUT)"),
        md("## Figure 1\n\nThe saved PNG and TIFF are the generated image used in the manuscript and submission files. Error bars describe uncertainty in the median, not variation among individuals."),
        code("make_main_figure(results['primary'], OUTPUT)\nplt.close('all')\ndisplay(Image(filename=str(OUTPUT / 'figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.png')))"),
        md("## Figure 2\n\nThe diagram separates changes in cohort membership from changes in measured trajectories. It contains no estimated effect sizes or causal decomposition."),
        code("make_conceptual_figure(OUTPUT)\nplt.close('all')\ndisplay(Image(filename=str(OUTPUT / 'figure2_conceptual_mechanism.png')))"),
        md("## Figure 1 source values"), code("show_outcomes(results['primary'])"),
        md("## Fixed-cohort sensitivity analysis\n\nThe following paired comparison reports the conditions quoted in the main manuscript. The primary analysis reselects participants at each detector setting; the sensitivity analysis retains the perfect-detector cohort. Neither analysis imposes the 25-day seizure-free interval rule."),
        code("show_main_cohort_comparison(results)"),
        md("## Validation"), code("validate_artifacts(results)")]
    return nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    notebook = build_notebook()
    if args.execute:
        NotebookClient(notebook, timeout=1200, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}}).execute()
    output = ROOT / "for_rtm3_paper.ipynb"
    nbf.write(notebook, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
