"""Tables, figures, and narrative summaries shared by the RTM3 notebooks."""

from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from IPython.display import Markdown, display


FILES = {
    "primary": "rtm3_primary_results.csv",
    "interval": "rtm3_interval_sensitivity_results.csv",
    "correction": "rtm3_correction_sensitivity_results.csv",
    "selection": "rtm3_selection_diagnostics.csv",
    "cohort": "rtm3_cohort_sensitivity_results.csv",
    "threshold": "rtm3_threshold_sensitivity_results.csv",
}


def load_results(directory=Path("results")):
    directory = Path(directory)
    results = {key: pd.read_csv(directory / filename) for key, filename in FILES.items()}
    for key, filename in [("metadata", "rtm3_run_metadata.json"), ("validation", "rtm3_validation.json")]:
        results[key] = json.loads((directory / filename).read_text())
    return results


def load_paired_results(directory=Path("results/paired_type3")):
    """Load the new replay results; notebook execution never runs the simulation."""
    directory = Path(directory)
    metadata = json.loads((directory / "paired_metadata.json").read_text())
    for filename, expected in metadata["artifact_sha256"].items():
        if hashlib.sha256((directory / filename).read_bytes()).hexdigest() != expected:
            raise AssertionError(f"Paired artifact hash mismatch: {filename}")
    return {"summary": pd.read_csv(directory / "paired_summary.csv"), "metadata": metadata,
            "validation": json.loads((directory / "paired_validation.json").read_text()),
            "directory": directory}


def paired_condition_labels(frame):
    return [f"s={s:.1f}; FAR={f:.1f}/day" for s, f in zip(frame.sensitivity, frame.FAR)]


def show_paired_results(paired):
    """Every requested condition, four paired cells, both denominators and diagnostics."""
    frame = paired["summary"]
    labels = paired_condition_labels(frame)
    display(Markdown("### Paired classification counts"))
    display(pd.DataFrame({"Condition": labels, "Selected n": frame.n_selected,
        "Both RTM n": frame.n_both_RTM, "Measurement-only n": frame.n_measurement_only,
        "Latent-only n": frame.n_latent_only, "Neither RTM n": frame.n_neither_RTM
    }).style.hide(axis="index"))
    display(Markdown("### Paired prevalence and discordance"))
    display(pd.DataFrame({"Condition": labels,
        "Measured RTM %": 100*frame.measured_RTM_fraction,
        "Latent RTM %": 100*frame.latent_RTM_fraction,
        "Latent-only / selected %": 100*frame.latent_only_fraction_selected,
        "Measured minus latent, pp": 100*frame.net_measured_minus_latent_fraction
    }).style.hide(axis="index").format(precision=2))
    display(Markdown("### Measurement-only fractions and conditional 95% Wilson intervals"))
    display(pd.DataFrame({"Condition": labels, "Measured RTM n": frame.n_measured_RTM,
        "Measurement-only / selected %, 95% CI": [f"{100*v:.2f} ({100*lo:.2f} to {100*hi:.2f})" for v,lo,hi in zip(frame.measurement_only_fraction_selected,frame.measurement_only_selected_ci_lo,frame.measurement_only_selected_ci_hi)],
        "Measurement-only / measured RTM %, 95% CI": [f"{100*v:.2f} ({100*lo:.2f} to {100*hi:.2f})" for v,lo,hi in zip(frame.measurement_only_fraction_measured_RTM,frame.measurement_only_measured_RTM_ci_lo,frame.measurement_only_measured_RTM_ci_hi)]
    }).style.hide(axis="index"))
    display(Markdown("### Retained latent eligibility and baseline diagnostics"))
    display(pd.DataFrame({"Condition": labels, "Selected n": frame.n_selected,
        "Latent eligible n": frame.n_selected_latent_eligible,
        "Latent ineligible n": frame.n_selected_latent_ineligible,
        "Zero true baseline n": frame.n_selected_zero_true_baseline
    }).style.hide(axis="index"))
    low, middle, perfect, alarm_low, alarm_high = [row for _, row in frame.iterrows()]
    display(Markdown(f"""At 10% sensitivity and no false alarms, {low.n_measurement_only:,.0f} of {low.n_selected:,.0f} selected participants had measurement-only RTM ({100*low.measurement_only_fraction_selected:.2f}%). They represented {100*low.measurement_only_fraction_measured_RTM:.2f}% of the {low.n_measured_RTM:,.0f} participants with measured RTM. Reverse discordance affected {low.n_latent_only:,.0f} participants ({100*low.latent_only_fraction_selected:.2f}%), leaving a net measured-minus-latent difference of {100*low.net_measured_minus_latent_fraction:.2f} percentage points. At 50% sensitivity, the primary measurement-only fraction was {100*middle.measurement_only_fraction_selected:.2f}% and the net difference was {100*middle.net_measured_minus_latent_fraction:.2f} points.

At FAR 0.1/day with perfect sensitivity, {alarm_low.n_measurement_only:,.0f} of {alarm_low.n_selected:,.0f} selected participants had measurement-only RTM ({100*alarm_low.measurement_only_fraction_selected:.2f}%). At FAR 1.0/day, {alarm_high.n_measurement_only:,.0f} of {alarm_high.n_selected:,.0f} did ({100*alarm_high.measurement_only_fraction_selected:.2f}%); these accounted for {100*alarm_high.measurement_only_fraction_measured_RTM:.2f}% of measured RTM. The latter condition also had {alarm_high.n_latent_only:,.0f} latent-only cases ({100*alarm_high.latent_only_fraction_selected:.2f}%), so its net prevalence difference was {100*alarm_high.net_measured_minus_latent_fraction:.2f} points. Perfect detection had {perfect.n_both_RTM:,.0f} RTM cases on both scales and no discordant cases.

The FAR 0.1/day and 1.0/day cohorts retained {alarm_low.n_selected_latent_ineligible:,.0f} and {alarm_high.n_selected_latent_ineligible:,.0f} participants who failed latent eligibility, including {alarm_low.n_selected_zero_true_baseline:,.0f} and {alarm_high.n_selected_zero_true_baseline:,.0f} with zero true baseline. These participants contribute to every full-selected denominator. Their inclusion is part of the detector-selected question, not a missing-data exception.

The gross measurement-only category and reverse discordance describe two directions of paired classification change. Their difference gives the net prevalence change within a condition; neither quantity estimates an exclusive biological mechanism. Comparisons between detector conditions also reflect changing enrollment. No MPC endpoint is calculated in this paired analysis."""))


def show_paired_reproducibility(paired):
    metadata = paired["metadata"]
    validation = paired["validation"]
    assert validation["status"] == "PASS"
    assert metadata["all_21_released_counts_match"]
    display(Markdown("The separate replay command generated the new participant audit. This notebook reads its saved results and checks their hashes; it does not repeat latent generation or detector replay. The original released metadata below remains unchanged. Full array hashes, detector-grid hashes, environment differences, and per-condition validation are retained in `results/paired_type3/paired_metadata.json` and its companion audit files.\n\n```bash\npython rtm3_paired.py\npython scripts/validate_paired_type3.py\npython -m unittest discover -s tests -p test_rtm3_paired.py -v\n```"))
    print(json.dumps({"paired_analysis": metadata["analysis"], "source_commit": metadata["source_commit"],
        "seeds": {k:metadata[k] for k in ["latent_seed","sensitivity_detection_seed","false_alarm_seed"]},
        "source_array_hashes": metadata["source_array_hashes"],
        "source_file_hashes": metadata["source_file_hashes"],
        "environment_mismatches_from_release": metadata["environment_mismatches_from_release"],
        "match_scope": metadata["match_scope"],
        "independent_integer_classification_crosscheck": metadata["independent_integer_classification_crosscheck"],
        "selected_records_independently_checked": validation["records_checked"],
        "artifact_sha256": metadata["artifact_sha256"], "validation_status": validation["status"]}, indent=2))


def condition_labels(frame):
    return np.where(frame.sweep.eq("sensitivity"),
                    frame.sensitivity.map(lambda x: f"Sensitivity {100*x:.0f}%"),
                    frame.FAR.map(lambda x: f"FAR {x:.1f}/day"))


def show_outcomes(frame, group=None):
    """Print every condition, with denominators and confidence intervals."""
    groups = frame.groupby(group, sort=False) if group else [(None, frame)]
    for label, subset in groups:
        if label:
            display(Markdown(f"### {label}"))
        table = pd.DataFrame({
            "Condition": condition_labels(subset),
            "Selected n": subset.n_eligible,
            "Selected %": 100 * subset.eligibility_fraction,
            "Observed RTM %": 100 * subset.observed_RTM_fraction,
            "MPC %": subset.MPC_median,
            "MPC 95% CI": [f"{lo:.2f} to {hi:.2f}" for lo, hi in zip(subset.MPC_ci_lo, subset.MPC_ci_hi)],
            "MPC denominator": subset.n_valid_MPC,
            "Zero baseline n": subset.n_zero_baseline,
        })
        display(table.style.hide(axis="index").format({
            "Selected %": "{:.2f}", "Observed RTM %": "{:.2f}", "MPC %": "{:.2f}"}))


def show_selection(results):
    frame = results["primary"]
    for sweep, subset in frame.groupby("sweep", sort=False):
        display(Markdown(f"### {sweep} sweep"))
        table = pd.DataFrame({
            "Condition": condition_labels(subset),
            "Selected n": subset.n_eligible,
            "Selected %": 100 * subset.eligibility_fraction,
            "Median sampled rate": subset.median_sampled_monthly_frequency_eligible,
            "Median effective long-term rate": subset.median_effective_long_term_frequency_eligible,
            "Mean baseline": subset.mean_baseline_frequency_eligible,
            "Mean test": subset.mean_test_frequency_eligible,
            "Mean baseline minus reference": subset.mean_baseline_excess_over_detector_long_term_eligible,
            "Mean test minus reference": subset.mean_test_excess_over_detector_long_term_eligible,
        })
        display(table.style.hide(axis="index").format(precision=2))


def show_main_cohort_comparison(results):
    """Show the paired conditions quoted in the main manuscript."""
    frame = results["cohort"]
    selected = frame[
        (frame.sweep.eq("FAR") & frame.FAR.isin([0, 1]))
        | (frame.sweep.eq("sensitivity") & frame.sensitivity.eq(.1))
    ]
    assert len(selected) == 6
    show_outcomes(selected, "cohort_strategy")
    fixed = selected[selected.cohort_strategy.eq("Fixed perfect-detector cohort")]
    reference = fixed[fixed.sweep.eq("FAR") & fixed.FAR.eq(0)].iloc[0]
    low = fixed[fixed.sweep.eq("sensitivity")].iloc[0]
    high = fixed[fixed.sweep.eq("FAR") & fixed.FAR.eq(1)].iloc[0]
    display(Markdown(f"""The same {reference.n_eligible:,.0f} participants selected under perfect detection were retained in every fixed-cohort condition. Their latent seizure histories and detector realizations were shared with the reselected-cohort analysis.

Observed RTM was {100*reference.observed_RTM_fraction:.1f}% with perfect detection, {100*low.observed_RTM_fraction:.1f}% at 10% sensitivity, and {100*high.observed_RTM_fraction:.1f}% at one false alarm/day. Thus, the increases at these endpoints in the primary analysis were absent when membership was held fixed. This comparison is not an additive causal decomposition of RTM Type 3.

MPC denotes the median percentage reduction from baseline; positive values indicate a reduction. Fixed-cohort MPC was {low.MPC_median:.1f}% at 10% sensitivity and {high.MPC_median:.1f}% at one false alarm/day. It excluded {low.n_zero_baseline:,.0f} and {high.n_zero_baseline:,.0f} participants with zero measured baseline, respectively. Its denominator therefore changes with detector performance even though cohort membership is fixed. The Appendix provides every condition, the full selection diagnostics, and the remaining sensitivity analyses."""))


def comparison_figure(frame, group, title, output=None):
    """Plot selection and both outcomes, without implying a causal decomposition."""
    sweeps = list(frame.sweep.drop_duplicates())
    fig, axes = plt.subplots(len(sweeps), 3, figsize=(13.5, 4 * len(sweeps)), squeeze=False)
    for row, sweep in enumerate(sweeps):
        xcol, scale = ("sensitivity", 100) if sweep == "sensitivity" else ("FAR", 1)
        for index, (label, part) in enumerate(frame[frame.sweep.eq(sweep)].groupby(group, sort=False)):
            part = part.sort_values(xcol)
            for axis, metric, yscale, ylabel in zip(axes[row],
                    ["eligibility_fraction", "observed_RTM_fraction", "MPC_median"],
                    [100, 100, 1], ["Selected participants (%)", "Observed RTM (%)", "Placebo MPC (%)"]):
                axis.plot(part[xcol] * scale, part[metric] * yscale,
                          marker=["o", "s", "^", "D", "v", "x"][index % 6],
                          linestyle=["-", "--", "-.", ":", "--", "-."][index % 6],
                          markersize=3, label=label)
                axis.set_xlabel("Sensitivity (%)" if sweep == "sensitivity" else "False alarm rate (alarms/day)")
                axis.set_ylabel(ylabel)
                axis.grid(alpha=.25)
                axis.spines[["top", "right"]].set_visible(False)
        # Presentation-only limit for this figure's false-alarm selection panel.
        if group == "cohort_strategy" and title == "Fixed versus reselected cohorts" and sweep == "FAR":
            axes[row, 0].set_ylim(0, 40)
            axes[row, 0].set_yticks(np.arange(0, 41, 10))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(title, fontsize=15)
    fig.legend(handles, labels, loc="lower center", ncol=1 if len(labels) > 2 else 2, fontsize=9)
    bottom = .08 + .032 * len(labels) / len(sweeps)
    fig.tight_layout(rect=(0, bottom, 1, .96))
    if output:
        fig.savefig(output, dpi=180, bbox_inches="tight")
    return fig


def interval_comparison(results):
    original = results["interval"].query("eligibility_variant == 'Original 25-day rule retained'")
    paired = results["primary"].merge(original, on=["sweep", "sensitivity", "FAR"], suffixes=("_main", "_interval"), validate="one_to_one")
    table = pd.DataFrame({"Condition": condition_labels(paired)})
    for metric, scale, label in [("eligibility_fraction", 100, "Selection difference, pp"),
                                 ("observed_RTM_fraction", 100, "Observed RTM difference, pp"),
                                 ("MPC_median", 1, "MPC difference, pp")]:
        table[label] = scale * (paired[metric + "_interval"] - paired[metric + "_main"])
    display(table.style.hide(axis="index").format(precision=3))
    display(pd.DataFrame({"Metric": table.columns[1:], "Maximum absolute difference, pp": table.iloc[:, 1:].abs().max().values}).style.hide(axis="index").format(precision=3))
    combined = pd.concat([results["primary"].assign(rule="Main: monthly counts only"), original.assign(rule="Monthly counts plus original 25-day rule")])
    return comparison_figure(combined, "rule", "Sensitivity analysis of the 25-day eligibility rule")


def threshold_comparison(results):
    frame = results["threshold"]
    assert frame.n_eligibility_mask_disagreements.eq(0).all()
    a = frame[frame.threshold_strategy.eq("Primary corrected-count eligibility")].set_index("FAR")
    b = frame[frame.threshold_strategy.eq("Raw counts with additive expected-alarm thresholds")].set_index("FAR")
    columns = ["n_eligible", "n_observed_RTM", "observed_RTM_fraction", "MPC_median", "MPC_ci_lo", "MPC_ci_hi", "n_valid_MPC"]
    pd.testing.assert_frame_equal(a[columns], b[columns], check_exact=True)
    display(a[["each_baseline_month_threshold_raw", "mean_baseline_threshold_raw", "baseline_total_threshold_raw", "n_eligibility_mask_disagreements"]].style.format(precision=1))
    show_outcomes(frame, "threshold_strategy")
    display(Markdown("Across all 11 false-alarm conditions, additive raw-count thresholds selected exactly the same participants as corrected-count eligibility. All observed-RTM counts, median changes, and confidence limits were identical. This changes the expression of the rule, not its statistical precision. Outcome counts remain corrected in both analyses."))


def show_design(results):
    metadata = results["metadata"]
    keys = ["n_patients_per_condition", "baseline_months", "test_months", "long_term_frequency_months", "base_seed", "sensitivity_detection_seed", "false_alarm_seed"]
    display(pd.DataFrame({"Parameter": keys, "Value": [str(metadata[key]) for key in keys]}).style.hide(axis="index"))
    display(pd.DataFrame(metadata["software"].items(), columns=["Software", "Version"]).style.hide(axis="index"))
    display(Markdown("Simulator provenance and shared-patient diagnostics:"))
    display(metadata["simulator_provenance"])
    display(metadata["cohort_generation_diagnostics"])


def show_fixed_membership_diagnostics(results):
    """Report retained identities that fail condition-specific enrollment."""
    fixed = results["cohort"].query("cohort_strategy == 'Fixed perfect-detector cohort'")
    table = pd.DataFrame({
        "Condition": condition_labels(fixed),
        "Retained n": fixed.n_eligible,
        "Retained meeting current criteria n": fixed.n_eligible - fixed.n_selected_failing_current_eligibility,
        "Retained failing current criteria n": fixed.n_selected_failing_current_eligibility,
        "Eligible in all potential participants n": fixed.n_currently_eligible,
    })
    display(table.style.hide(axis="index"))


def show_recorded_metadata(results):
    """Preserve the released run record separately from the reporting build."""
    display(Markdown("```json\n" + json.dumps(results["metadata"], indent=2, sort_keys=True) + "\n```"))


def validate_artifacts(results, root=Path.cwd()):
    root = Path(root)
    metadata = results["metadata"]
    assert results["validation"]["status"] == "PASS"
    assert hashlib.sha256((root / "reproduce_rtm3.py").read_bytes()).hexdigest() == metadata["script_sha256"]
    assert hashlib.sha256((root / "realSim.py").read_bytes()).hexdigest() == metadata["simulator_provenance"]["local_source_sha256"]
    assert len(results["primary"]) == 21
    assert results["primary"].n_total.eq(metadata["n_patients_per_condition"]).all()
    display(results["validation"])
    display(Markdown("PASS: the current analysis and simulator files match the recorded source hashes. All 21 primary conditions use the recorded cohort size."))


def summarize_results(results):
    primary = results["primary"]
    perfect = primary.query("sweep == 'FAR' and FAR == 0").iloc[0]
    low = primary.query("sweep == 'sensitivity' and sensitivity == 0.1").iloc[0]
    high = primary.query("sweep == 'FAR' and FAR == 1").iloc[0]
    fixed = results["cohort"].query("cohort_strategy == 'Fixed perfect-detector cohort'")
    flo = fixed.query("sweep == 'sensitivity' and sensitivity == 0.1").iloc[0]
    fhi = fixed.query("sweep == 'FAR' and FAR == 1").iloc[0]
    raw = results["correction"].query("correction_strategy == 'Uncorrected observed counts; 25-day rule omitted' and FAR == 1").iloc[0]
    validity = results["validation"]
    display(Markdown(f"""With perfect detection, {perfect.n_eligible:,.0f} of {perfect.n_total:,.0f} participants were selected, observed RTM was {100*perfect.observed_RTM_fraction:.1f}%, and placebo median percentage change was {perfect.MPC_median:.1f}%. Under reselection, the corresponding observed-RTM and median-change estimates were {100*low.observed_RTM_fraction:.1f}% and {low.MPC_median:.1f}% at 10% sensitivity, and {100*high.observed_RTM_fraction:.1f}% and {high.MPC_median:.1f}% at one false alarm/day.

Holding the perfect-detector cohort fixed gave observed RTM of {100*flo.observed_RTM_fraction:.1f}% and median change of {flo.MPC_median:.1f}% at 10% sensitivity, and {100*fhi.observed_RTM_fraction:.1f}% and {fhi.MPC_median:.1f}% at one false alarm/day. Median change was defined for {flo.n_valid_MPC:,.0f} and {fhi.n_valid_MPC:,.0f} selected participants, respectively; {flo.n_zero_baseline:,.0f} and {fhi.n_zero_baseline:,.0f} had zero measured baseline. These comparisons describe the influence of allowing cohort membership to change. They are not additive causal components of RTM, and the fixed-cohort median has a performance-dependent valid-baseline denominator.

Moving expected false alarms from the counts to the eligibility thresholds produced exactly identical results. It does not reduce the remaining random false-alarm error.

Retaining the original 25-day rule changed selection by at most {validity['max_absolute_eligibility_difference_primary_vs_original_25day_pp']:.2f} percentage points, observed RTM by {validity['max_absolute_observed_RTM_difference_primary_vs_original_25day_pp']:.2f} points, and MPC by {validity['max_absolute_MPC_difference_primary_vs_original_25day_pp']:.2f} points. The main directions remained. Alternative interval rules gave different results because observed alarms and oracle seizure labels define different zero-day runs. The original rule's larger overall maxima include the sensitivity sweep; the observed-day and oracle-day alternatives were evaluated only on the 11-point FAR sweep. On that common FAR sweep, the observed-day and original rules have equal maximum absolute changes from the main analysis: 0.9571218755 percentage points for observed RTM and 0.8658008658 percentage points for MPC. The oracle-day rule has larger maxima on that same FAR grid. These comparisons concern maxima across conditions, not a uniform ordering at every condition.

Without any count correction or threshold adjustment, FAR one/day selected {raw.n_eligible:,.0f} participants, observed RTM was {100*raw.observed_RTM_fraction:.1f}%, and MPC was {raw.MPC_median:.1f}%. This does not establish that uncorrected counts are preferable: alarms allow low-seizure-rate participants to pass fixed raw thresholds and contribute to both the baseline denominator and test counts. The reference now includes expected alarms so that the RTM comparison stays on the observation scale. Exact oracle removal retains the perfect-detector results at every FAR because it removes all realized simulated alarms. It is a diagnostic, not an available correction for unlabeled clinical events."""))


def make_conceptual_figure(output_dir=None):
    """Draw the manuscript's conceptual diagram; no patient values are plotted."""
    fig, ax = plt.subplots(figsize=(12, 7.2))
    ax.set(xlim=(0, 12), ylim=(0, 7.2))
    ax.axis("off")
    ax.text(6, 6.95, "How detector error can change apparent placebo response", ha="center", va="top", fontsize=18, weight="bold")
    ax.text(6, 6.48, "Conceptual pathways, not a causal decomposition of the simulation", ha="center", fontsize=11, color="#475569")
    boxes = [
        (.25, 4.7, 3.3, 1.2, "Latent seizure history", "Variable baseline and test counts", "#e2e8f0"),
        (4.35, 4.7, 3.3, 1.2, "Detector observations", "Missed seizures and false alarms", "#dbeafe"),
        (8.45, 4.7, 3.3, 1.2, "Corrected monthly counts", "Expected alarms subtracted;\nnegative counts set to zero", "#dbeafe"),
        (2.25, 2.55, 3.6, 1.2, "Eligibility selection", "Baseline count thresholds\ncan change cohort membership", "#fef3c7"),
        (6.5, 2.55, 3.6, 1.2, "Measured trial trajectory", "Baseline and test counts\ncan change within a fixed cohort", "#fef3c7"),
        (3.75, .45, 4.5, 1.15, "Observed trial outcomes", "Observed-RTM indicator and\nmedian percentage change", "#dcfce7"),
    ]
    for x, y, width, height, title, body, color in boxes:
        ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.04", facecolor=color, edgecolor="#64748b", linewidth=1.2))
        ax.text(x + width/2, y + height-.25, title, ha="center", va="center", fontsize=12, weight="bold")
        ax.text(x + width/2, y + .42, body, ha="center", va="center", fontsize=10.5, linespacing=1.5)
    arrows = [((3.6, 5.3), (4.25, 5.3)), ((7.7, 5.3), (8.35, 5.3)),
              ((9.9, 4.6), (4.4, 3.85)), ((10.2, 4.6), (8.5, 3.85)),
              ((4.1, 2.45), (5.4, 1.7)), ((8.15, 2.45), (6.6, 1.7))]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#334155"))
    fig.subplots_adjust(0, 0, 1, 1)
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / "figure2_conceptual_mechanism.png", dpi=300, facecolor="white")
        fig.savefig(output_dir / "figure2_conceptual_mechanism.tif", dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    return fig
