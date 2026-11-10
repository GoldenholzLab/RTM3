#!/usr/bin/env python3
"""Build the methodological appendix using shared reporting functions."""
import argparse
import os
import sys
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import AsyncKernelManager

ROOT = Path(__file__).resolve().parents[1]


def execute_in_process(notebook):
    """Execute every code cell in order when the host forbids socket listeners.

    Uses a fresh IPython shell with the inline backend and captures rich display
    outputs in nbformat form. This is reporting execution, not a simulation run.
    """
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    from matplotlib_inline.backend_inline import configure_inline_support
    import matplotlib
    import matplotlib.pyplot as plt

    InteractiveShell.clear_instance()
    shell = InteractiveShell.instance()
    matplotlib.use("module://matplotlib_inline.backend_inline")
    configure_inline_support(shell, "inline")
    previous_directory = Path.cwd()
    previous_path = list(sys.path)
    try:
        os.chdir(ROOT)
        sys.path.insert(0, str(ROOT))
        for cell in notebook.cells:
            if cell.cell_type != "code":
                continue
            cell.outputs = []
            with capture_output() as captured:
                result = shell.run_cell(cell.source, store_history=True)
            if result.error_before_exec or result.error_in_exec:
                raise result.error_before_exec or result.error_in_exec
            cell.execution_count = result.execution_count
            if captured.stdout:
                cell.outputs.append(nbf.v4.new_output("stream", name="stdout", text=captured.stdout))
            if captured.stderr:
                cell.outputs.append(nbf.v4.new_output("stream", name="stderr", text=captured.stderr))
            for output in captured.outputs:
                cell.outputs.append(nbf.v4.new_output("display_data", data=output.data, metadata=output.metadata))
        notebook.metadata["reporting_execution"] = {
            "engine": "IPython InteractiveShell in process",
            "order": "All code cells, top to bottom in a fresh shell",
            "simulation_rerun": False,
        }
    finally:
        plt.close("all")
        os.chdir(previous_directory)
        sys.path[:] = previous_path
        InteractiveShell.clear_instance()
    return notebook


def build_notebook():
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    cells = [md((ROOT / "appendix_methods.md").read_text()),
        code("from pathlib import Path\nimport matplotlib.pyplot as plt\nfrom rtm3_reporting import *\nresults = load_results(Path('results'))\nshow_design(results)"),
        md("## Primary results\n\nThe sensitivity sweep changes sensitivity from 10% to 100% with no false alarms. The false-alarm sweep changes alarms from 0 to 1 per day in steps of 0.1 at 100% sensitivity. These are two separate sweeps, not a full factorial grid. Their perfect-detector reference is repeated, giving 21 rows for 20 distinct detector settings. Latent histories, treatment effect (zero), monthly correction, and eligibility thresholds stay constant; eligibility is reapplied to each condition. Selection percentages use all 100,000 potential participants. Observed RTM uses the selected denominator, and MPC uses selected participants with positive baseline. Every condition is reported. Positive MPC indicates apparent improvement without treatment."),
        code("show_outcomes(results['primary'])"),
        md("## Selection diagnostics\n\nThis analysis asks whether changing detector performance changes who enters the primary cohort and how baseline and test counts differ from the reference. Frequencies are seizures per 30-day month among participants meeting primary criteria. Sampled and effective long-term rates are latent; baseline, test, and excess-over-reference measures are on the corrected detector scale. Median effective long-term frequency among selected participants is 8.42 under perfect detection, 39.24 at 10% sensitivity, and 7.39 at one false alarm/day. Selection therefore changes the underlying frequency distribution. At these settings, mean baseline excess over reference is 1.32, 2.73, and 3.60 per month, while mean test excess is -0.08, -0.11, and 0.60. These selected-cohort summaries describe associations and do not identify causal contributions."),
        code("show_selection(results)"),
        md("## Fixed and reselected cohorts\n\n**The fixed perfect-detector cohort is hypothetical and completely unrealistic as an enrollment strategy for a trial whose eligibility is determined by the imperfect detector being evaluated.**\n\nSome retained participants would have failed that detector's entry criteria. The comparison tests how the result depends on changing enrollment. Reselection lets the current detector determine entry; fixed membership retains the 35,393 perfect-detector identities. Latent histories and the exact detector realizations are shared. The table of membership diagnostics distinguishes retained patients who meet or fail the current criteria from eligibility among all 100,000 potential patients. MPC applies only to positive-baseline participants; zero baselines are reported separately.\n\nPerfect detection gives observed RTM of 37.2% and MPC of 14.3% in both strategies. At 10% sensitivity, reselection gives 64.4% observed RTM and 48.1% MPC, compared with 29.3% and 33.3% under fixed membership. At one false alarm/day, the corresponding values are 53.6% and 31.2% under reselection, versus 34.7% and 10.4% under fixed membership. The reversal relative to perfect detection concerns observed-RTM prevalence: it rises under reselection and falls under fixed membership at these endpoints. MPC does not reverse in every condition; fixed-cohort MPC at 10% sensitivity still exceeds 14.3%. The comparison is hypothetical, is not an additive causal decomposition, and does not isolate RTM Type 3."),
        code("show_outcomes(results['cohort'], 'cohort_strategy')\ncomparison_figure(results['cohort'], 'cohort_strategy', 'Fixed versus reselected cohorts')\nplt.show()\nshow_fixed_membership_diagnostics(results)"),
        md("## Additive expected-alarm thresholds\n\nThis check asks whether adding expected alarms to the raw eligibility thresholds differs from subtracting them from the monthly counts. Only eligibility uses raw counts; outcomes, references, latent histories, and detector realizations retain the primary definitions. Every monthly minimum and mean threshold must move together. The equality proof above and zero mask disagreements at all 11 false-alarm rates establish equivalence. This re-expression does not remove random alarm noise or increase precision."),
        code("threshold_comparison(results)"),
        md("## The 25-day rule\n\nThis analysis asks whether the primary findings depend on omitting a maximum 25-day seizure-free interval. It changes only the interval requirement while preserving corrected monthly thresholds, outcomes, detector realizations, and references. Differences below are sensitivity analysis minus main analysis, in percentage points. The original rounded-day rule changes selection by at most 2.284 points, observed RTM by 2.492 points, and MPC by 2.288 points. Its discontinuity at 0.6 alarms/day comes from daily rounding. The alternatives use observed days or oracle detected-seizure days. The released label \"25-day rule on latent true days\" refers to the oracle daily-count alternative in the false-alarm sweep, where sensitivity is 100% and detected latent seizures equal latent seizures. Those two alternatives are evaluated at the 11 false-alarm settings only, whereas the original rule is evaluated across both sweeps."),
        code("show_outcomes(results['interval'], 'eligibility_variant')\ninterval_comparison(results)\nplt.show()"),
        md("## Correction and interval alternatives\n\nThis comparison asks how results depend on the count scale and on how seizure-free intervals are defined. Sensitivity stays at 100%; false alarms vary from 0 to 1 per day. All strategies share latent histories and condition-specific detector events. Expected monthly correction is the primary method. The original, observed-day, and oracle-day interval alternatives retain corrected monthly counts and change only the 25-day criterion. Uncorrected outcomes use raw counts for both enrollment and outcomes, keep the original fixed raw thresholds, and add 30 times the daily false-alarm rate to the long-term reference. Exact oracle removal subtracts every realized labeled false alarm and uses the latent reference. It is unavailable for unlabeled clinical events. None of these alternatives establishes which real detector or correction is preferable. At one false alarm/day, the main, original-interval, and observed-interval strategies each select 35,134 participants and yield 53.6% observed RTM and 31.2% MPC. The oracle-day interval rule selects 27,003 and yields 49.0% and 26.7%. Raw counts select all 100,000 and yield 27.1% and 0.0%, whereas exact oracle removal retains 35,393 and yields 37.2% and 14.3%. These endpoints illustrate why the interval definition and observation scale must accompany each result."),
        code("show_outcomes(results['correction'], 'correction_strategy')\ncomparison_figure(results['correction'], 'correction_strategy', 'False-alarm correction and interval alternatives')\nplt.show()"),
        md((ROOT / "paired_type3_methods.md").read_text()),
        code("paired = load_paired_results(Path('results/paired_type3'))\nshow_paired_results(paired)"),
        md("## Interpretation"), code("summarize_results(results)"),
        md("## Participant-level replay and audit provenance"), code("show_paired_reproducibility(paired)"),
        md("## Reproducibility checks\n\nThe following checks read released aggregates and verify the recorded simulation source hashes. They do not rerun the simulation. The recorded validation status and all reported validation results are retained."), code("validate_artifacts(results)"),
        md("## Recorded simulation settings\n\nThis is the complete released run metadata. Its timestamp, software versions, source path, and artifact hashes describe the original simulation, not this remote presentation build. The separate reporting audit records the environment used to execute this notebook and generate the appendix."),
        code("show_recorded_metadata(results)")]
    return nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--kernel-transport", choices=["tcp", "ipc", "inprocess"], default="tcp",
                        help="Use inprocess when the host forbids all socket listeners.")
    args = parser.parse_args()
    notebook = build_notebook()
    if args.execute and args.kernel_transport == "inprocess":
        execute_in_process(notebook)
    elif args.execute:
        manager = AsyncKernelManager(kernel_name="python3", transport=args.kernel_transport)
        client = NotebookClient(notebook, km=manager, timeout=1200, resources={"metadata": {"path": str(ROOT)}})
        try:
            client.execute()
        finally:
            client._cleanup_kernel()
    nbf.validate(notebook)
    output = ROOT / "for_appendix_RTM3.ipynb"
    nbf.write(notebook, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
