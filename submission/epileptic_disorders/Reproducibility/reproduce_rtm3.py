#!/usr/bin/env python3
"""Canonical, deterministic reproduction of the RTM3 simulation analyses.

The primary analysis uses only the two monthly seizure-count eligibility
criteria.  The historical 25-day seizure-free-interval criterion is evaluated
as a sensitivity analysis, including alternative ways to handle false alarms
for a day-level rule.

Run from the repository root:

    python reproduce_rtm3.py

The script writes all numerical results, diagnostics, metadata, and the main
manuscript figure to ``results/`` by default.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats

from realSim import get_mSF, simulator_base


DEFAULT_SEED = 20260813
DEFAULT_N_PATIENTS = 100_000

DAYS_PER_MONTH = 30
BASELINE_MONTHS = 2
TEST_MONTHS = 3
TOTAL_MONTHS = BASELINE_MONTHS + TEST_MONTHS
TOTAL_DAYS = TOTAL_MONTHS * DAYS_PER_MONTH
LONG_TERM_MONTHS = 36

ELIGIBILITY_MEAN_MIN = 4.0
ELIGIBILITY_MONTHLY_MIN = 3.0
ELIGIBILITY_MAX_ZERO_RUN = 25

SENSITIVITY_VALUES = np.round(np.arange(0.1, 1.01, 0.1), 1)
FAR_VALUES = np.round(np.arange(0.0, 1.01, 0.1), 1)

SIMULATOR_PROVENANCE = {
    "upstream_repository": "https://github.com/GoldenholzLab/CHOCOLATE",
    "upstream_commit": "03fdf198d176028e1a495e109029e45f980ba67e",
    "upstream_source_path": "/Users/dgoldenh/Documents/GitHub/CHOCOLATE/realSim.py",
    "source_file": "realSim.py",
    "previous_source_sha256": "7c4f1ac426a136b7539ed6ab868164e3af67d2322cb9ed558fe34f5428e3ee9d",
    "upstream_source_sha256": "7c4f1ac426a136b7539ed6ab868164e3af67d2322cb9ed558fe34f5428e3ee9d",
    "source_migration": "The requested upstream file is byte-identical to the previous local file.",
}

RESELECTED_COHORT = "Reselected at each detector condition"
FIXED_COHORT = "Fixed perfect-detector cohort"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-patients",
        type=int,
        default=DEFAULT_N_PATIENTS,
        help=f"Potential participants per condition (default: {DEFAULT_N_PATIENTS:,}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Base random seed (default: {DEFAULT_SEED}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV, JSON, PNG, and TIFF outputs (default: results).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10_000,
        help="Progress interval during latent-cohort generation; 0 disables it.",
    )
    return parser.parse_args()


def monthly_counts(daily_counts: np.ndarray) -> np.ndarray:
    """Aggregate a patient-by-day array into 30-day months."""
    daily_counts = np.asarray(daily_counts)
    if daily_counts.ndim != 2 or daily_counts.shape[1] != TOTAL_DAYS:
        raise ValueError(f"Expected shape (n, {TOTAL_DAYS}), got {daily_counts.shape}")
    return daily_counts.reshape(len(daily_counts), TOTAL_MONTHS, DAYS_PER_MONTH).sum(axis=2)


def longest_zero_run(daily_counts: np.ndarray) -> np.ndarray:
    """Longest consecutive zero-count run for every patient."""
    daily_counts = np.asarray(daily_counts)
    current = np.zeros(len(daily_counts), dtype=np.int16)
    longest = np.zeros(len(daily_counts), dtype=np.int16)
    for day_index in range(daily_counts.shape[1]):
        is_zero = daily_counts[:, day_index] == 0
        current = np.where(is_zero, current + 1, 0)
        longest = np.maximum(longest, current)
    return longest


def eligible_mask(
    monthly: np.ndarray,
    daily_for_interval: np.ndarray | None = None,
    include_25_day_rule: bool = False,
) -> np.ndarray:
    """Apply the primary monthly criteria and, optionally, the 25-day rule."""
    baseline_monthly = np.asarray(monthly)[:, :BASELINE_MONTHS]
    eligible = (
        (baseline_monthly.mean(axis=1) >= ELIGIBILITY_MEAN_MIN)
        & (baseline_monthly.min(axis=1) >= ELIGIBILITY_MONTHLY_MIN)
    )
    if include_25_day_rule:
        if daily_for_interval is None:
            raise ValueError("Daily counts are required when the 25-day rule is included.")
        baseline_daily = np.asarray(daily_for_interval)[:, : BASELINE_MONTHS * DAYS_PER_MONTH]
        eligible &= longest_zero_run(baseline_daily) <= ELIGIBILITY_MAX_ZERO_RUN
    return eligible


def additive_threshold_mask(observed_monthly: np.ndarray, far: float) -> np.ndarray:
    """Select on raw counts after adding expected alarms to each monthly threshold."""
    baseline = np.asarray(observed_monthly)[:, :BASELINE_MONTHS]
    expected_monthly_alarms = DAYS_PER_MONTH * far
    return (
        (baseline.mean(axis=1) >= ELIGIBILITY_MEAN_MIN + expected_monthly_alarms)
        & (baseline.min(axis=1) >= ELIGIBILITY_MONTHLY_MIN + expected_monthly_alarms)
    )


def observed_rtm_indicator(
    monthly: np.ndarray,
    effective_monthly_rate: np.ndarray,
    sensitivity: float,
    reference_offset_monthly: float = 0.0,
    long_term_total: np.ndarray | None = None,
) -> np.ndarray:
    """Evaluate the strict endpoint with exact scaled integer counts.

    The canonical grid has integer monthly counts and rational sensitivities.
    The reference is an integer 36-month seizure total divided by 36. Retaining
    that numerator avoids classifying equal baseline and test distances as RTM.
    Noninteger counts outside this grid use float64 comparisons with an explicit
    1e-10 relative tolerance for equality.
    """
    monthly = np.asarray(monthly, dtype=np.float64)
    effective = np.asarray(effective_monthly_rate, dtype=np.float64)
    if long_term_total is None:
        # Backward compatibility for callers supplying only a 36-month mean.
        long_term_total = np.rint(effective * LONG_TERM_MONTHS).astype(np.int64)
    else:
        long_term_total = np.asarray(long_term_total, dtype=np.int64)
    sensitivity_fraction = Fraction(str(float(sensitivity))).limit_denominator(10_000)
    offset_fraction = Fraction(str(float(reference_offset_monthly))).limit_denominator(10_000)
    reference_denominator = LONG_TERM_MONTHS * sensitivity_fraction.denominator
    scale = math.lcm(BASELINE_MONTHS, TEST_MONTHS, reference_denominator, offset_fraction.denominator)
    integer_counts = np.all(np.abs(monthly - np.rint(monthly)) <= 1e-10)
    represents_reference = np.allclose(
        effective, long_term_total / LONG_TERM_MONTHS, rtol=1e-7, atol=1e-9
    )
    if integer_counts and represents_reference and scale <= 10**9:
        counts = np.rint(monthly).astype(np.int64)
        baseline = counts[:, :BASELINE_MONTHS].sum(axis=1) * (scale // BASELINE_MONTHS)
        test = counts[:, BASELINE_MONTHS:].sum(axis=1) * (scale // TEST_MONTHS)
        reference = (
            long_term_total * sensitivity_fraction.numerator * (scale // reference_denominator)
            + offset_fraction.numerator * (scale // offset_fraction.denominator)
        )
        return (baseline > reference) & (baseline - reference > np.abs(test - reference))
    baseline = monthly[:, :BASELINE_MONTHS].mean(axis=1)
    test = monthly[:, BASELINE_MONTHS:].mean(axis=1)
    reference = effective * sensitivity + reference_offset_monthly
    tolerance = 1e-10 * np.maximum.reduce([np.ones_like(reference), np.abs(baseline), np.abs(test), np.abs(reference)])
    return (baseline > reference + tolerance) & (
        baseline - reference > np.abs(test - reference) + tolerance
    )


def wilson_interval(successes: int, total: int, confidence: float = 0.95) -> tuple[float, float]:
    if total == 0:
        return np.nan, np.nan
    z_value = stats.norm.ppf(0.5 + confidence / 2)
    estimate = successes / total
    denominator = 1 + z_value**2 / total
    center = (estimate + z_value**2 / (2 * total)) / denominator
    half_width = z_value * np.sqrt(
        estimate * (1 - estimate) / total + z_value**2 / (4 * total**2)
    ) / denominator
    return float(center - half_width), float(center + half_width)


def exact_median_interval(
    values: np.ndarray, confidence: float = 0.95
) -> tuple[float, float, float]:
    values = np.sort(np.asarray(values, dtype=float))
    values = values[np.isfinite(values)]
    n_values = len(values)
    if n_values == 0:
        return np.nan, np.nan, np.nan
    alpha = 1 - confidence
    lower_index = int(stats.binom.ppf(alpha / 2, n_values, 0.5)) - 1
    upper_index = int(stats.binom.isf(alpha / 2, n_values, 0.5))
    lower = -np.inf if lower_index < 0 else float(values[lower_index])
    upper = np.inf if upper_index >= n_values else float(values[upper_index])
    return float(np.median(values)), lower, upper


def summarize_condition(
    *,
    monthly: np.ndarray,
    sampled_monthly_rate: np.ndarray,
    effective_monthly_rate: np.ndarray,
    sensitivity: float,
    far: float,
    sweep: str,
    eligibility_variant: str,
    daily_for_interval: np.ndarray | None = None,
    include_25_day_rule: bool = False,
    selection_mask: np.ndarray | None = None,
    cohort_strategy: str = RESELECTED_COHORT,
    reference_offset_monthly: float = 0.0,
    long_term_total: np.ndarray | None = None,
) -> dict[str, float | int | str | bool]:
    """Summarize outcomes and selection diagnostics for one condition."""
    monthly = np.asarray(monthly, dtype=float)
    currently_eligible = eligible_mask(
        monthly,
        daily_for_interval=daily_for_interval,
        include_25_day_rule=include_25_day_rule,
    )
    eligible = currently_eligible if selection_mask is None else np.asarray(selection_mask, dtype=bool)
    if eligible.shape != currently_eligible.shape:
        raise ValueError("The selection mask must contain one entry per latent participant.")
    n_total = len(eligible)
    n_eligible = int(eligible.sum())
    baseline = monthly[:, :BASELINE_MONTHS].mean(axis=1)
    test = monthly[:, BASELINE_MONTHS:].mean(axis=1)

    detector_specific_long_term = np.asarray(effective_monthly_rate, dtype=np.float64) * sensitivity + reference_offset_monthly
    observed_rtm = observed_rtm_indicator(
        monthly, effective_monthly_rate, sensitivity, reference_offset_monthly, long_term_total
    )
    valid_mpc = eligible & (baseline > 0) & np.isfinite(baseline) & np.isfinite(test)
    n_valid_mpc = int(valid_mpc.sum())
    n_zero_baseline = int((eligible & (baseline == 0)).sum())

    if n_eligible:
        percentage_change = 100 * (1 - test[valid_mpc] / baseline[valid_mpc])
        mpc, mpc_lo, mpc_hi = exact_median_interval(percentage_change)
        n_observed_rtm = int(observed_rtm[eligible].sum())
        observed_rtm_fraction = n_observed_rtm / n_eligible
        rtm_lo, rtm_hi = wilson_interval(n_observed_rtm, n_eligible)
    else:
        mpc = mpc_lo = mpc_hi = np.nan
        n_observed_rtm = 0
        observed_rtm_fraction = rtm_lo = rtm_hi = np.nan

    eligibility_lo, eligibility_hi = wilson_interval(n_eligible, n_total)

    def selected_mean(values: np.ndarray) -> float:
        return float(np.mean(values[eligible])) if n_eligible else np.nan

    def selected_median(values: np.ndarray) -> float:
        return float(np.median(values[eligible])) if n_eligible else np.nan

    return {
        "sweep": sweep,
        "sensitivity": float(sensitivity),
        "FAR": float(far),
        "eligibility_variant": eligibility_variant,
        "cohort_strategy": cohort_strategy,
        "includes_25_day_rule": bool(include_25_day_rule),
        "n_total": n_total,
        "n_eligible": n_eligible,
        "n_currently_eligible": int(currently_eligible.sum()),
        "n_selected_failing_current_eligibility": int((eligible & ~currently_eligible).sum()),
        "n_valid_MPC": n_valid_mpc,
        "n_zero_baseline": n_zero_baseline,
        "n_undefined_MPC": n_eligible - n_valid_mpc,
        "MPC_valid_fraction_of_selected": n_valid_mpc / n_eligible if n_eligible else np.nan,
        "reference_offset_monthly": reference_offset_monthly,
        "eligibility_fraction": n_eligible / n_total,
        "eligibility_ci_lo": eligibility_lo,
        "eligibility_ci_hi": eligibility_hi,
        "n_observed_RTM": n_observed_rtm,
        "observed_RTM_fraction": observed_rtm_fraction,
        "observed_RTM_ci_lo": rtm_lo,
        "observed_RTM_ci_hi": rtm_hi,
        "MPC_median": mpc,
        "MPC_ci_lo": mpc_lo,
        "MPC_ci_hi": mpc_hi,
        "median_sampled_monthly_frequency_eligible": selected_median(
            np.asarray(sampled_monthly_rate)
        ),
        "median_effective_long_term_frequency_eligible": selected_median(
            np.asarray(effective_monthly_rate)
        ),
        "median_detector_specific_long_term_frequency_eligible": selected_median(
            detector_specific_long_term
        ),
        "mean_baseline_frequency_eligible": selected_mean(baseline),
        "mean_test_frequency_eligible": selected_mean(test),
        "median_baseline_frequency_eligible": selected_median(baseline),
        "median_test_frequency_eligible": selected_median(test),
        "mean_baseline_excess_over_detector_long_term_eligible": selected_mean(
            baseline - detector_specific_long_term
        ),
        "mean_test_excess_over_detector_long_term_eligible": selected_mean(
            test - detector_specific_long_term
        ),
    }


def generate_patient_histories(sampled_rate: float) -> tuple[tuple, tuple]:
    """Draw separate histories with the same participant rate, cycles, and cluster status."""
    trial = simulator_base(
        sampRATE=1, number_of_days=TOTAL_DAYS,
        defaultSeizureFreq=sampled_rate, returnDetails=True,
    )
    frequencies, amplitudes, cluster_status = trial[5], trial[6], bool(trial[8])
    reference = simulator_base(
        sampRATE=1, number_of_days=LONG_TERM_MONTHS * DAYS_PER_MONTH,
        defaultSeizureFreq=sampled_rate,
        cyclesTF=len(frequencies) > 0,
        CP=[np.array(frequencies, copy=True), np.array(amplitudes, copy=True)],
        clusterParams=[float(cluster_status), 1, 7, 1, 1],
        returnDetails=True,
    )
    if bool(reference[8]) != cluster_status:
        raise AssertionError("The reference changed participant cluster status.")
    if not np.array_equal(frequencies, reference[5]) or not np.array_equal(amplitudes, reference[6]):
        raise AssertionError("The reference changed participant cycle parameters.")
    return trial, reference


def generate_latent_cohort(
    n_patients: int, seed: int, progress_every: int, diagnostics: dict | None = None,
    long_term_totals_out: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Generate a common latent cohort and a separate 36-month frequency estimate."""
    np.random.seed(seed)  # realSim uses NumPy's legacy global random state.
    latent_daily = np.empty((n_patients, TOTAL_DAYS), dtype=np.int16)
    sampled_monthly_rate = np.empty(n_patients, dtype=np.float64)
    effective_monthly_rate = np.empty(n_patients, dtype=np.float64)
    long_term_totals = np.empty(n_patients, dtype=np.int64) if long_term_totals_out is None else long_term_totals_out
    if long_term_totals.shape != (n_patients,) or long_term_totals.dtype != np.int64:
        raise ValueError("long_term_totals_out must be an int64 array with one entry per participant.")
    cluster_status = np.empty(n_patients, dtype=bool)
    nonzero_cycles = np.empty(n_patients, dtype=bool)

    started = time.perf_counter()
    for patient_index in range(n_patients):
        sampled_rate = float(np.asarray(get_mSF(requested_msf=-1)).reshape(-1)[0])
        sampled_monthly_rate[patient_index] = sampled_rate
        trial, reference = generate_patient_histories(sampled_rate)
        if np.any(trial[0] < 0) or np.any(trial[0] > np.iinfo(np.int16).max):
            raise AssertionError("Latent counts exceed the storage range.")
        latent_daily[patient_index] = trial[0].astype(np.int16)
        long_term_totals[patient_index] = int(reference[0].sum())
        effective_monthly_rate[patient_index] = long_term_totals[patient_index] / LONG_TERM_MONTHS
        cluster_status[patient_index] = bool(trial[8])
        nonzero_cycles[patient_index] = np.any(trial[7])
        completed = patient_index + 1
        if progress_every and completed % progress_every == 0:
            print(f"Generated {completed:,}/{n_patients:,} latent participants", flush=True)

    if diagnostics is not None:
        diagnostics.update({
            "shared_rate_cluster_status_and_cycle_parameters": True,
            "n_clustered_participants": int(cluster_status.sum()),
            "n_participants_with_nonzero_daily_cycles": int(nonzero_cycles.sum()),
            "long_term_total_sha256": hashlib.sha256(long_term_totals.tobytes()).hexdigest(),
            "cycle_phase_handling": "Independent phases for separate histories; frequency and amplitude arrays are fixed per participant.",
        })
    return (
        latent_daily,
        sampled_monthly_rate,
        effective_monthly_rate,
        time.perf_counter() - started,
    )


def run_sensitivity_sweep(
    *,
    latent_daily: np.ndarray,
    sampled_monthly_rate: np.ndarray,
    effective_monthly_rate: np.ndarray,
    seed: int,
    long_term_total: np.ndarray | None = None,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Nested sensitivity sweep using the same latent patients in every condition."""
    detector_rng = np.random.default_rng(seed)
    detected_daily = np.zeros_like(latent_daily, dtype=np.int16)
    undetected_daily = latent_daily.copy()
    previous_sensitivity = 0.0
    primary_rows: list[dict] = []
    sensitivity_rows: list[dict] = []
    cohort_rows: list[dict] = []
    fixed_selection = eligible_mask(monthly_counts(latent_daily))
    if long_term_total is None:
        long_term_total = np.rint(effective_monthly_rate * LONG_TERM_MONTHS).astype(np.int64)

    for sensitivity in SENSITIVITY_VALUES:
        conditional_probability = (
            (float(sensitivity) - previous_sensitivity) / (1 - previous_sensitivity)
            if previous_sensitivity < 1
            else 0.0
        )
        newly_detected = detector_rng.binomial(
            undetected_daily, conditional_probability
        ).astype(np.int16)
        detected_daily += newly_detected
        undetected_daily -= newly_detected
        detected_monthly = monthly_counts(detected_daily)

        primary = summarize_condition(
                monthly=detected_monthly,
                sampled_monthly_rate=sampled_monthly_rate,
                effective_monthly_rate=effective_monthly_rate,
                sensitivity=float(sensitivity),
                far=0.0,
                sweep="sensitivity",
                eligibility_variant="Primary monthly criteria",
                long_term_total=long_term_total,
            )
        primary_rows.append(primary)
        cohort_rows.append(primary)
        cohort_rows.append(summarize_condition(
            monthly=detected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=float(sensitivity), far=0.0, sweep="sensitivity",
            eligibility_variant="Perfect-detector monthly criteria held fixed",
            selection_mask=fixed_selection, cohort_strategy=FIXED_COHORT,
            long_term_total=long_term_total,
        ))
        sensitivity_rows.append(
            summarize_condition(
                monthly=detected_monthly,
                sampled_monthly_rate=sampled_monthly_rate,
                effective_monthly_rate=effective_monthly_rate,
                sensitivity=float(sensitivity),
                far=0.0,
                sweep="sensitivity",
                eligibility_variant="Original 25-day rule retained",
                daily_for_interval=detected_daily,
                include_25_day_rule=True,
                long_term_total=long_term_total,
            )
        )
        previous_sensitivity = float(sensitivity)

    if np.any(undetected_daily):
        raise AssertionError("The 100% sensitivity condition did not detect all latent seizures.")
    return primary_rows, sensitivity_rows, cohort_rows


def run_far_sweep(
    *,
    latent_daily: np.ndarray,
    sampled_monthly_rate: np.ndarray,
    effective_monthly_rate: np.ndarray,
    seed: int,
    long_term_total: np.ndarray | None = None,
) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    """Nested FAR sweep, primary monthly correction, and interval sensitivities."""
    false_alarm_rng = np.random.default_rng(seed)
    cumulative_false_alarms = np.zeros_like(latent_daily, dtype=np.int16)
    latent_monthly = monthly_counts(latent_daily)
    previous_far = 0.0
    primary_rows: list[dict] = []
    interval_rows: list[dict] = []
    correction_rows: list[dict] = []
    cohort_rows: list[dict] = []
    threshold_rows: list[dict] = []
    fixed_selection = eligible_mask(latent_monthly)
    if long_term_total is None:
        long_term_total = np.rint(effective_monthly_rate * LONG_TERM_MONTHS).astype(np.int64)

    for far in FAR_VALUES:
        far = float(far)
        if far > 0:
            far_increment = far - previous_far
            cumulative_false_alarms += false_alarm_rng.poisson(
                far_increment, size=latent_daily.shape
            ).astype(np.int16)
        observed_daily = latent_daily + cumulative_false_alarms
        observed_monthly = monthly_counts(observed_daily)
        expected_monthly_alarms = far * DAYS_PER_MONTH
        corrected_monthly = np.maximum(observed_monthly - expected_monthly_alarms, 0.0)

        primary = summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="Primary monthly criteria",
            long_term_total=long_term_total,
        )
        primary_rows.append(primary)
        cohort_rows.append(primary)
        cohort_rows.append(summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0, far=far, sweep="FAR",
            eligibility_variant="Perfect-detector monthly criteria held fixed",
            selection_mask=fixed_selection, cohort_strategy=FIXED_COHORT,
            long_term_total=long_term_total,
        ))
        raw_threshold_selection = additive_threshold_mask(observed_monthly, far)
        mask_disagreements = int(np.count_nonzero(raw_threshold_selection != eligible_mask(corrected_monthly)))
        threshold_common = {
            "mean_baseline_threshold_raw": ELIGIBILITY_MEAN_MIN + expected_monthly_alarms,
            "each_baseline_month_threshold_raw": ELIGIBILITY_MONTHLY_MIN + expected_monthly_alarms,
            "baseline_total_threshold_raw": BASELINE_MONTHS * (ELIGIBILITY_MEAN_MIN + expected_monthly_alarms),
            "n_eligibility_mask_disagreements": mask_disagreements,
            "outcome_count_basis": "Expected monthly correction with zero floor",
        }
        threshold_rows.append({
            **primary, **threshold_common,
            "threshold_strategy": "Primary corrected-count eligibility",
        })
        raw_threshold_summary = summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0, far=far, sweep="FAR",
            eligibility_variant="Raw counts with additive expected-alarm thresholds",
            selection_mask=raw_threshold_selection, long_term_total=long_term_total,
        )
        threshold_rows.append({
            **raw_threshold_summary, **threshold_common,
            "threshold_strategy": "Raw counts with additive expected-alarm thresholds",
        })
        correction_rows.append(
            {
                **primary,
                "correction_strategy": "Expected monthly correction; 25-day rule omitted",
            }
        )

        daily_subtraction = int(np.rint(far))
        original_daily_corrected = np.maximum(observed_daily - daily_subtraction, 0)
        original_interval = summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="Original 25-day rule retained",
            daily_for_interval=original_daily_corrected,
            include_25_day_rule=True,
            long_term_total=long_term_total,
        )
        interval_rows.append(original_interval)
        correction_rows.append(
            {
                **original_interval,
                "correction_strategy": "Expected monthly correction + original daily rounding",
            }
        )

        observed_interval = summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="25-day rule on observed days",
            daily_for_interval=observed_daily,
            include_25_day_rule=True,
            long_term_total=long_term_total,
        )
        interval_rows.append(observed_interval)
        correction_rows.append(
            {
                **observed_interval,
                "correction_strategy": "Expected monthly correction + observed daily intervals",
            }
        )

        oracle_interval = summarize_condition(
            monthly=corrected_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="25-day rule on latent true days",
            daily_for_interval=latent_daily,
            include_25_day_rule=True,
            long_term_total=long_term_total,
        )
        interval_rows.append(oracle_interval)
        correction_rows.append(
            {
                **oracle_interval,
                "correction_strategy": "Expected monthly correction + oracle daily intervals",
            }
        )

        uncorrected = summarize_condition(
            monthly=observed_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="Monthly criteria without FAR correction",
            reference_offset_monthly=expected_monthly_alarms,
            long_term_total=long_term_total,
        )
        correction_rows.append(
            {
                **uncorrected,
                "correction_strategy": "Uncorrected observed counts; 25-day rule omitted",
            }
        )

        exact_oracle = summarize_condition(
            monthly=latent_monthly,
            sampled_monthly_rate=sampled_monthly_rate,
            effective_monthly_rate=effective_monthly_rate,
            sensitivity=1.0,
            far=far,
            sweep="FAR",
            eligibility_variant="Exact oracle false-alarm removal",
            long_term_total=long_term_total,
        )
        correction_rows.append(
            {
                **exact_oracle,
                "correction_strategy": "Exact oracle false-alarm removal",
            }
        )
        previous_far = far

    return primary_rows, interval_rows, correction_rows, cohort_rows, threshold_rows


def make_main_figure(primary_results: pd.DataFrame, output_dir: Path | None = None) -> plt.Figure:
    sensitivity = primary_results[primary_results.sweep == "sensitivity"].sort_values(
        "sensitivity"
    )
    far = primary_results[primary_results.sweep == "FAR"].sort_values("FAR")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 18,
            "axes.linewidth": 1.4,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "lines.linewidth": 2.5,
        }
    )
    color = "#1976B9"
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex="col")

    axes[0, 0].plot(
        100 * sensitivity.sensitivity,
        100 * sensitivity.observed_RTM_fraction,
        marker="o",
        color=color,
    )
    axes[0, 1].plot(
        far.FAR,
        100 * far.observed_RTM_fraction,
        marker="o",
        color=color,
    )

    sensitivity_yerr = np.vstack(
        [
            sensitivity.MPC_median - sensitivity.MPC_ci_lo,
            sensitivity.MPC_ci_hi - sensitivity.MPC_median,
        ]
    )
    far_yerr = np.vstack(
        [far.MPC_median - far.MPC_ci_lo, far.MPC_ci_hi - far.MPC_median]
    )
    axes[1, 0].errorbar(
        100 * sensitivity.sensitivity,
        sensitivity.MPC_median,
        yerr=sensitivity_yerr,
        marker="o",
        capsize=4,
        color=color,
    )
    axes[1, 1].errorbar(
        far.FAR,
        far.MPC_median,
        yerr=far_yerr,
        marker="o",
        capsize=4,
        color=color,
    )

    perfect_mpc = float(
        sensitivity.loc[np.isclose(sensitivity.sensitivity, 1.0), "MPC_median"].iloc[0]
    )
    for axis in axes[1]:
        axis.axhline(perfect_mpc, color="#777777", linestyle="--", linewidth=1.5)

    for axis in axes.flat:
        axis.grid(True, color="#D9D9D9", linestyle="--", linewidth=0.8, alpha=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[0, 0].set_ylabel("Observed RTM (%)")
    axes[1, 0].set_ylabel("Placebo MPC (%)")
    axes[1, 0].set_xlabel("Sensitivity (%)")
    axes[1, 1].set_xlabel("False alarm rate (alarms/day)")
    axes[0, 0].set_ylim(0, 100)
    axes[0, 1].set_ylim(0, 100)
    axes[1, 0].set_ylim(bottom=min(-5, axes[1, 0].get_ylim()[0]))
    axes[1, 1].set_ylim(bottom=min(-5, axes[1, 1].get_ylim()[0]))

    axes[0, 0].set_title("A", fontweight="bold")
    axes[0, 1].set_title("B", fontweight="bold")
    axes[1, 0].set_title("C", fontweight="bold")
    axes[1, 1].set_title("D", fontweight="bold")
    fig.suptitle(
        "Detector performance changes observed RTM and apparent placebo response",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.012,
        "Primary eligibility uses monthly seizure-count criteria only; dashed line is the perfect-detector MPC.",
        ha="center",
        fontsize=11,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.09, top=0.90, hspace=0.34, wspace=0.16)

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / "figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.png", dpi=300, facecolor="white")
        fig.savefig(
            output_dir / "figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.tif",
            dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"},
        )
    return fig


def make_selection_figure(primary_results: pd.DataFrame, output_dir: Path) -> None:
    sensitivity = primary_results[primary_results.sweep == "sensitivity"].sort_values(
        "sensitivity"
    )
    far = primary_results[primary_results.sweep == "FAR"].sort_values("FAR")
    color = "#1976B9"
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex="col")
    axes[0, 0].plot(
        100 * sensitivity.sensitivity,
        100 * sensitivity.eligibility_fraction,
        marker="o",
        color=color,
    )
    axes[0, 1].plot(
        far.FAR, 100 * far.eligibility_fraction, marker="o", color=color
    )
    axes[1, 0].plot(
        100 * sensitivity.sensitivity,
        sensitivity.median_effective_long_term_frequency_eligible,
        marker="o",
        color=color,
    )
    axes[1, 1].plot(
        far.FAR,
        far.median_effective_long_term_frequency_eligible,
        marker="o",
        color=color,
    )
    axes[0, 0].set_ylabel("Eligible (%)")
    axes[1, 0].set_ylabel("Median effective long-term\nfrequency among eligible")
    axes[1, 0].set_xlabel("Sensitivity (%)")
    axes[1, 1].set_xlabel("False alarm rate (alarms/day)")
    axes[0, 0].set_title("Sensitivity sweep")
    axes[0, 1].set_title("False-alarm sweep")
    for axis in axes.flat:
        axis.grid(True, color="#D9D9D9", linestyle="--", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle("Primary-analysis selection diagnostics", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_dir / "appendix_selection_diagnostics.png", dpi=200, facecolor="white")
    plt.close(fig)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_results(
    primary_results: pd.DataFrame,
    original_interval_results: pd.DataFrame,
    correction_results: pd.DataFrame,
    cohort_results: pd.DataFrame | None = None,
    threshold_results: pd.DataFrame | None = None,
) -> dict[str, float | str | bool]:
    sensitivity_perfect = primary_results[
        (primary_results.sweep == "sensitivity")
        & np.isclose(primary_results.sensitivity, 1.0)
    ].iloc[0]
    far_zero = primary_results[
        (primary_results.sweep == "FAR") & np.isclose(primary_results.FAR, 0.0)
    ].iloc[0]
    comparison_columns = [
        "n_eligible",
        "observed_RTM_fraction",
        "MPC_median",
        "MPC_ci_lo",
        "MPC_ci_hi",
    ]
    perfect_condition_match = all(
        np.isclose(sensitivity_perfect[column], far_zero[column])
        for column in comparison_columns
    )

    original = original_interval_results[
        original_interval_results.eligibility_variant == "Original 25-day rule retained"
    ]
    paired = primary_results.merge(
        original[
            ["sweep", "sensitivity", "FAR", "observed_RTM_fraction", "MPC_median", "eligibility_fraction"]
        ],
        on=["sweep", "sensitivity", "FAR"],
        suffixes=("_primary", "_25day"),
        validate="one_to_one",
    )
    max_mpc_difference = float(
        np.max(np.abs(paired.MPC_median_primary - paired.MPC_median_25day))
    )
    max_rtm_difference_pp = float(
        100
        * np.max(
            np.abs(
                paired.observed_RTM_fraction_primary - paired.observed_RTM_fraction_25day
            )
        )
    )
    max_eligibility_difference_pp = float(
        100
        * np.max(
            np.abs(paired.eligibility_fraction_primary - paired.eligibility_fraction_25day)
        )
    )

    oracle = correction_results[
        correction_results.correction_strategy == "Exact oracle false-alarm removal"
    ]
    oracle_invariant = bool(
        oracle[["n_eligible", "observed_RTM_fraction", "MPC_median"]].nunique().eq(1).all()
    )
    mpc_denominators_complete = bool(
        (primary_results.n_valid_MPC + primary_results.n_undefined_MPC).eq(primary_results.n_eligible).all()
    )
    fixed_cohort_size_invariant = True
    perfect_fixed_match = True
    if cohort_results is not None:
        fixed = cohort_results[cohort_results.cohort_strategy == FIXED_COHORT]
        fixed_cohort_size_invariant = bool(fixed.n_eligible.eq(int(far_zero.n_eligible)).all())
        mpc_denominators_complete &= bool(
            (cohort_results.n_valid_MPC + cohort_results.n_undefined_MPC).eq(cohort_results.n_eligible).all()
            and cohort_results.n_undefined_MPC.eq(cohort_results.n_zero_baseline).all()
        )
        fixed_perfect = fixed[(fixed.sensitivity == 1.0) & (fixed.FAR == 0.0)]
        perfect_fixed_match = all(
            np.isclose(fixed_perfect[column], far_zero[column]).all() for column in comparison_columns
        )
    additive_thresholds_equivalent = True
    if threshold_results is not None:
        additive_thresholds_equivalent = bool(threshold_results.n_eligibility_mask_disagreements.eq(0).all())
        for _, rows in threshold_results.groupby("FAR"):
            additive_thresholds_equivalent &= bool(rows[comparison_columns].nunique(dropna=False).eq(1).all())

    return {
        "status": "PASS"
        if (perfect_condition_match and oracle_invariant and primary_results.notna().all().all()
            and mpc_denominators_complete and fixed_cohort_size_invariant
            and perfect_fixed_match and additive_thresholds_equivalent)
        else "FAIL",
        "perfect_detector_condition_identical_between_sweeps": bool(perfect_condition_match),
        "exact_oracle_invariant_across_FAR": oracle_invariant,
        "primary_rows": int(len(primary_results)),
        "fixed_cohort_size_invariant": fixed_cohort_size_invariant,
        "fixed_and_reselected_perfect_detector_match": perfect_fixed_match,
        "MPC_denominators_complete": mpc_denominators_complete,
        "additive_raw_thresholds_equivalent_to_primary": additive_thresholds_equivalent,
        "max_absolute_MPC_difference_primary_vs_original_25day_pp": max_mpc_difference,
        "max_absolute_observed_RTM_difference_primary_vs_original_25day_pp": max_rtm_difference_pp,
        "max_absolute_eligibility_difference_primary_vs_original_25day_pp": max_eligibility_difference_pp,
    }


def run_analysis(
    *, n_patients: int, seed: int, output_dir: Path, progress_every: int
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"RTM3 canonical run: N={n_patients:,}; base seed={seed}")

    cohort_generation_diagnostics: dict = {}
    long_term_total = np.empty(n_patients, dtype=np.int64)
    latent_daily, sampled_rate, effective_rate, latent_seconds = generate_latent_cohort(
        n_patients=n_patients,
        seed=seed,
        progress_every=progress_every,
        diagnostics=cohort_generation_diagnostics,
        long_term_totals_out=long_term_total,
    )
    print(f"Latent cohort generated in {latent_seconds:.1f} seconds")

    sensitivity_primary, sensitivity_interval, sensitivity_cohort = run_sensitivity_sweep(
        latent_daily=latent_daily,
        sampled_monthly_rate=sampled_rate,
        effective_monthly_rate=effective_rate,
        seed=seed + 1,
        long_term_total=long_term_total,
    )
    far_primary, far_interval, correction_rows, far_cohort, threshold_rows = run_far_sweep(
        latent_daily=latent_daily,
        sampled_monthly_rate=sampled_rate,
        effective_monthly_rate=effective_rate,
        seed=seed + 2,
        long_term_total=long_term_total,
    )

    primary_results = pd.DataFrame(sensitivity_primary + far_primary)
    interval_results = pd.DataFrame(sensitivity_interval + far_interval)
    correction_results = pd.DataFrame(correction_rows)
    cohort_results = pd.DataFrame(sensitivity_cohort + far_cohort)
    threshold_results = pd.DataFrame(threshold_rows)
    diagnostics_columns = [
        "sweep",
        "sensitivity",
        "FAR",
        "n_total",
        "n_eligible",
        "n_valid_MPC",
        "n_zero_baseline",
        "n_undefined_MPC",
        "eligibility_fraction",
        "median_sampled_monthly_frequency_eligible",
        "median_effective_long_term_frequency_eligible",
        "median_detector_specific_long_term_frequency_eligible",
        "mean_baseline_frequency_eligible",
        "mean_test_frequency_eligible",
        "median_baseline_frequency_eligible",
        "median_test_frequency_eligible",
        "mean_baseline_excess_over_detector_long_term_eligible",
        "mean_test_excess_over_detector_long_term_eligible",
        "observed_RTM_fraction",
        "MPC_median",
    ]
    selection_diagnostics = primary_results[diagnostics_columns].copy()

    paths = {
        "primary": output_dir / "rtm3_primary_results.csv",
        "interval": output_dir / "rtm3_interval_sensitivity_results.csv",
        "correction": output_dir / "rtm3_correction_sensitivity_results.csv",
        "selection": output_dir / "rtm3_selection_diagnostics.csv",
        "cohort": output_dir / "rtm3_cohort_sensitivity_results.csv",
        "threshold": output_dir / "rtm3_threshold_sensitivity_results.csv",
        "validation": output_dir / "rtm3_validation.json",
        "metadata": output_dir / "rtm3_run_metadata.json",
        "figure_png": output_dir / "figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.png",
        "figure_tif": output_dir / "figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.tif",
        "selection_figure": output_dir / "appendix_selection_diagnostics.png",
    }
    primary_results.to_csv(paths["primary"], index=False)
    interval_results.to_csv(paths["interval"], index=False)
    correction_results.to_csv(paths["correction"], index=False)
    selection_diagnostics.to_csv(paths["selection"], index=False)
    cohort_results.to_csv(paths["cohort"], index=False)
    threshold_results.to_csv(paths["threshold"], index=False)
    main_figure = make_main_figure(primary_results, output_dir)
    plt.close(main_figure)
    make_selection_figure(primary_results, output_dir)

    validation = validate_results(primary_results, interval_results, correction_results, cohort_results, threshold_results)
    write_json(paths["validation"], validation)
    if validation["status"] != "PASS":
        raise AssertionError(f"Validation failed: {validation}")

    artifact_hashes = {
        key: sha256_file(path)
        for key, path in paths.items()
        if key not in {"metadata"} and path.exists()
    }
    metadata = {
        "analysis": "RTM3 canonical reproduction",
        "analysis_revision": "2026-09-15 shared participant parameters and exact RTM endpoint",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "primary_eligibility": {
            "mean_baseline_monthly_minimum": ELIGIBILITY_MEAN_MIN,
            "each_baseline_month_minimum": ELIGIBILITY_MONTHLY_MIN,
            "seizure_free_interval_rule": None,
        },
        "sensitivity_analysis_eligibility": {
            "maximum_seizure_free_interval_days": ELIGIBILITY_MAX_ZERO_RUN
        },
        "n_patients_per_condition": n_patients,
        "base_seed": seed,
        "latent_seed": seed,
        "sensitivity_detection_seed": seed + 1,
        "false_alarm_seed": seed + 2,
        "baseline_months": BASELINE_MONTHS,
        "test_months": TEST_MONTHS,
        "long_term_frequency_months": LONG_TERM_MONTHS,
        "long_term_reference": (
            "Independent 36-month history with the trial participant's sampled rate, cluster status, "
            "and cycle frequency/amplitude arrays. Integer total retained for exact RTM comparisons."
        ),
        "reference_count_scale": {
            "primary": "Sensitivity times true long-term frequency is the biological target; it excludes zero-floor bias.",
            "uncorrected_comparator": "Sensitivity times true long-term frequency plus 30 times FAR.",
        },
        "RTM_ties": "Exact scaled integer comparisons on the canonical grid; equal distances do not satisfy the strict endpoint.",
        "fixed_cohort_analysis": {
            "selection": "Perfect-detector primary monthly eligibility, reused at every condition.",
            "RTM_denominator": "All selected participants, including zero measured baseline counts.",
            "MPC_denominator": "Selected participants with positive measured baseline, reported as n_valid_MPC.",
            "undefined_MPC": "Zero measured baseline counts are reported as n_zero_baseline and n_undefined_MPC.",
        },
        "additive_threshold_analysis": {
            "selection": "Raw monthly mean >= 4 + 30*FAR and each raw baseline month >= 3 + 30*FAR.",
            "equivalent_total_baseline_threshold": "8 + 60*FAR over the two baseline months.",
            "outcomes": "Primary expected-monthly-corrected counts with the zero floor.",
        },
        "cohort_generation_diagnostics": cohort_generation_diagnostics,
        "simulator_provenance": {
            **SIMULATOR_PROVENANCE,
            "local_source_sha256": sha256_file(Path(__file__).resolve().with_name("realSim.py")),
        },
        "sensitivity_values": [float(value) for value in SENSITIVITY_VALUES],
        "FAR_values_per_day": [float(value) for value in FAR_VALUES],
        "FAR_correction": "Subtract expected false alarms per 30-day month, then floor at zero.",
        "common_random_numbers": (
            "All device conditions use the same latent cohort. Sensitivity and FAR realizations "
            "are nested within their respective sweeps."
        ),
        "runtime_seconds_latent_generation": latent_seconds,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "artifact_sha256": artifact_hashes,
        "validation": validation,
    }
    write_json(paths["metadata"], metadata)

    print("Validation: PASS")
    print(json.dumps(validation, indent=2))
    print(f"Outputs written to {output_dir.resolve()}")
    return paths


def main() -> int:
    args = parse_args()
    if args.n_patients <= 0:
        raise SystemExit("--n-patients must be positive")
    run_analysis(
        n_patients=args.n_patients,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_every=args.progress_every,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
