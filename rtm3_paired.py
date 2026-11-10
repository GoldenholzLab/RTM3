"""Paired measured/latent strict RTM within detector-selected primary cohorts.

Replay the released nested detector grids, then retain five requested conditions.
This module never calls the canonical run_analysis or writes its output files.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import reproduce_rtm3 as canonical

ROOT = Path(__file__).resolve().parent
EXPECTED_LONG_TERM_SHA = "2a66b529678d29fd1de585ee3bbb015c6d89ac44fd097aca3b2852fc896e0d5f"
CONDITIONS = [(1, 0), (5, 0), (10, 0), (10, 1), (10, 10)]
EXPECTED_SELECTED = [2180, 17865, 35393, 35453, 35134]


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_record(array):
    return {"shape": list(array.shape), "dtype": str(array.dtype),
            "sha256": hashlib.sha256(array.tobytes()).hexdigest()}


def software_record():
    return {"python": platform.python_version(), "platform": platform.platform(),
            **{p: importlib.metadata.version(p) for p in ["numpy", "scipy", "pandas", "matplotlib"]}}


def strict_from_totals(baseline_total, test_total, long_term_total,
                       sensitivity_numerator=1, sensitivity_denominator=1):
    """Independent exact integer cross-check using two strict test bounds.

    B = baseline_total/2; T = test_total/3; reference = a*L/(36*b).
    B > reference and |T-reference| < B-reference are equivalent to
    2*test_total < 3*baseline_total and
    9*b*baseline_total + 6*b*test_total > a*L.
    These two strict bounds also imply baseline above reference.
    Equality at any boundary is false. Inputs are nonnegative integer totals.
    """
    arrays = [np.asarray(x) for x in (baseline_total, test_total, long_term_total)]
    if any(not np.issubdtype(x.dtype, np.integer) or np.any(x < 0) for x in arrays):
        raise ValueError("Nonnegative integer totals are required.")
    a, b = int(sensitivity_numerator), int(sensitivity_denominator)
    if not (0 <= a <= b and b > 0):
        raise ValueError("Sensitivity must be a rational number from zero to one.")
    baseline, test, reference_total = [x.astype(np.int64) for x in arrays]
    if any(np.any(x > np.iinfo(np.int64).max // (100 * b)) for x in arrays):
        raise ValueError("Counts exceed safe int64 cross-product range.")
    return ((2*test < 3*baseline)
            & (9*b*baseline + 6*b*test > a*reference_total))


def paired_counts(measured_rtm, latent_rtm):
    measured, latent = np.asarray(measured_rtm, dtype=bool), np.asarray(latent_rtm, dtype=bool)
    if measured.shape != latent.shape or measured.ndim != 1:
        raise ValueError("Paired indicators must be one-dimensional and have equal length.")
    both = int(np.count_nonzero(measured & latent))
    measured_only = int(np.count_nonzero(measured & ~latent))
    latent_only = int(np.count_nonzero(~measured & latent))
    neither = int(np.count_nonzero(~measured & ~latent))
    n = len(measured)
    measured_n, latent_n = both + measured_only, both + latent_only
    primary_ci = canonical.wilson_interval(measured_only, n)
    secondary_ci = canonical.wilson_interval(measured_only, measured_n)
    ratio = lambda numerator, denominator: numerator / denominator if denominator else float("nan")
    return {"n_selected": n, "n_both_RTM": both, "n_measurement_only": measured_only,
            "n_latent_only": latent_only, "n_neither_RTM": neither,
            "n_measured_RTM": measured_n, "n_latent_RTM": latent_n,
            "measured_RTM_fraction": ratio(measured_n, n),
            "latent_RTM_fraction": ratio(latent_n, n),
            "measurement_only_fraction_selected": ratio(measured_only, n),
            "measurement_only_fraction_measured_RTM": ratio(measured_only, measured_n),
            "measurement_only_selected_ci_lo": primary_ci[0],
            "measurement_only_selected_ci_hi": primary_ci[1],
            "measurement_only_measured_RTM_ci_lo": secondary_ci[0],
            "measurement_only_measured_RTM_ci_hi": secondary_ci[1],
            "latent_only_fraction_selected": ratio(latent_only, n),
            "net_measured_minus_latent_fraction": ratio(measured_only-latent_only, n)}


def classify_condition(true_monthly, measured_monthly, long_term_total, sensitivity):
    """Eligibility is measured primary monthly eligibility; no true-count filter."""
    selected = canonical.eligible_mask(measured_monthly, include_25_day_rule=False)
    effective = long_term_total / canonical.LONG_TERM_MONTHS
    measured = canonical.observed_rtm_indicator(measured_monthly, effective, sensitivity,
                                               long_term_total=long_term_total)
    latent = canonical.observed_rtm_indicator(true_monthly, effective, 1.0,
                                             long_term_total=long_term_total)
    tenth = round(sensitivity * 10)
    independent_measured = strict_from_totals(measured_monthly[:, :2].sum(1).astype(np.int64),
        measured_monthly[:, 2:].sum(1).astype(np.int64), long_term_total, tenth, 10)
    independent_latent = strict_from_totals(true_monthly[:, :2].sum(1),
        true_monthly[:, 2:].sum(1), long_term_total)
    if not np.array_equal(measured, independent_measured) or not np.array_equal(latent, independent_latent):
        raise AssertionError("Canonical and independent integer classification disagree.")
    return selected, measured, latent


def load_or_generate_latent(cache, seed, n_patients):
    cache = Path(cache)
    metadata_path = cache.with_name("latent_replay_metadata.json")
    if cache.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata["seed"] != seed or metadata["n_patients"] != n_patients:
            raise AssertionError("Cache seed or population mismatch.")
        arrays = dict(np.load(cache, allow_pickle=False))
        for name, array in arrays.items():
            if array_record(array) != metadata["arrays"][name]:
                raise AssertionError(f"Cached source array changed: {name}")
        for name, sha in metadata["source_sha256"].items():
            if file_hash(ROOT / name) != sha:
                raise AssertionError(f"Generation source changed: {name}")
        return arrays, metadata
    total = np.empty(n_patients, dtype=np.int64)
    diagnostics = {}
    daily, sampled, effective, seconds = canonical.generate_latent_cohort(
        n_patients, seed, 10000, diagnostics, total)
    arrays = {"latent_daily": daily, "sampled_monthly_rate": sampled,
              "effective_monthly_rate": effective, "long_term_total": total}
    metadata = {"seed": seed, "n_patients": n_patients, "seconds": seconds,
                "diagnostics": diagnostics, "software": software_record(),
                "arrays": {k: array_record(a) for k, a in arrays.items()},
                "source_sha256": {f: file_hash(ROOT/f) for f in ["realSim.py", "reproduce_rtm3.py"]}}
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **arrays)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return arrays, metadata


def replay_detector_conditions(latent_daily, seed):
    """Yield all 21 original rows in original RNG order, including both controls."""
    rng = np.random.default_rng(seed+1)
    detected = np.zeros_like(latent_daily, dtype=np.int16)
    undetected = latent_daily.copy()
    previous = 0.0
    for sensitivity in canonical.SENSITIVITY_VALUES:
        probability = (float(sensitivity)-previous)/(1-previous) if previous < 1 else 0.0
        newly = rng.binomial(undetected, probability).astype(np.int16)
        detected += newly
        undetected -= newly
        monthly = canonical.monthly_counts(detected)
        yield "sensitivity", float(sensitivity), 0.0, detected, monthly, monthly
        previous = float(sensitivity)
    if np.any(undetected):
        raise AssertionError("Perfect sensitivity failed to recover all latent events.")
    rng = np.random.default_rng(seed+2)
    alarms = np.zeros_like(latent_daily, dtype=np.int16)
    previous = 0.0
    for far in canonical.FAR_VALUES:
        far = float(far)
        if far > 0:
            alarms += rng.poisson(far-previous, size=latent_daily.shape).astype(np.int16)
        observed_daily = latent_daily + alarms
        observed = canonical.monthly_counts(observed_daily)
        measured = np.maximum(observed - far*canonical.DAYS_PER_MONTH, 0.0)
        yield "FAR", 1.0, far, observed_daily, observed, measured
        previous = far


def run_paired_analysis(output_dir, cache):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    released_paths = list((ROOT/"results").glob("*.csv")) + list((ROOT/"results").glob("*.json"))
    released_hashes = {str(p.relative_to(ROOT)): file_hash(p) for p in released_paths}
    release = json.loads((ROOT/"results/rtm3_run_metadata.json").read_text())
    released = pd.read_csv(ROOT/"results/rtm3_primary_results.csv")
    arrays, latent_metadata = load_or_generate_latent(cache, release["base_seed"], release["n_patients_per_condition"])
    total, daily = arrays["long_term_total"], arrays["latent_daily"]
    if array_record(total)["sha256"] != EXPECTED_LONG_TERM_SHA:
        raise AssertionError("Long-term source hash does not match release; replay is not identical.")
    true = canonical.monthly_counts(daily)
    rows, records, grid_checks, detector_hashes = {}, {}, [], []
    for sweep, sensitivity, far, observed_daily, observed, measured in replay_detector_conditions(daily, release["base_seed"]):
        selected, measured_rtm, latent_rtm = classify_condition(true, measured, total, sensitivity)
        released_row = released[(released.sweep == sweep) & (released.sensitivity == sensitivity) & (released.FAR == far)].iloc[0]
        check = {"sweep": sweep, "sensitivity": sensitivity, "FAR": far,
                 "n_selected": int(selected.sum()), "n_measured_RTM": int(measured_rtm[selected].sum()),
                 "expected_selected": int(released_row.n_eligible), "expected_measured_RTM": int(released_row.n_observed_RTM)}
        check["exact_match"] = check["n_selected"] == check["expected_selected"] and check["n_measured_RTM"] == check["expected_measured_RTM"]
        grid_checks.append(check)
        (output/"replay_grid_validation.json").write_text(json.dumps(grid_checks, indent=2)+"\n")
        if not check["exact_match"]:
            raise AssertionError(f"Replay differs from canonical counts: {check}")
        detector_hashes.append({"sweep": sweep, "sensitivity": sensitivity, "FAR": far,
                               "observed_daily": array_record(observed_daily), "corrected_monthly": array_record(measured),
                               "selected_ids": array_record(np.flatnonzero(selected))})
        key = (round(sensitivity*10), round(far*10))
        if key in CONDITIONS:
            summary = {"condition_id": CONDITIONS.index(key), "sensitivity": sensitivity, "FAR": far,
                       **paired_counts(measured_rtm[selected], latent_rtm[selected])}
            latent_eligible = canonical.eligible_mask(true)
            summary.update({"n_selected_latent_eligible": int(np.count_nonzero(selected & latent_eligible)),
                            "n_selected_latent_ineligible": int(np.count_nonzero(selected & ~latent_eligible)),
                            "n_selected_zero_true_baseline": int(np.count_nonzero(selected & (true[:,:2].sum(1) == 0)))})
            if key in rows:
                if summary != rows[key]:
                    raise AssertionError("Shared perfect control differs between sweeps.")
            else:
                rows[key] = summary
                records[key] = {"participant_id": np.flatnonzero(selected).astype(np.int32),
                    "condition_id": np.full(selected.sum(), CONDITIONS.index(key), dtype=np.int8),
                    "sensitivity_tenths": np.full(selected.sum(), key[0], dtype=np.int8),
                    "FAR_tenths": np.full(selected.sum(), key[1], dtype=np.int8),
                    "true_monthly": true[selected].astype(np.int32),
                    "observed_monthly": observed[selected].astype(np.int32),
                    "measured_monthly": np.rint(measured[selected]).astype(np.int32),
                    "true_baseline_total": true[selected,:2].sum(1), "true_test_total": true[selected,2:].sum(1),
                    "measured_baseline_total": np.rint(measured[selected,:2].sum(1)).astype(np.int64),
                    "measured_test_total": np.rint(measured[selected,2:].sum(1)).astype(np.int64),
                    "long_term_total": total[selected], "measured_RTM": measured_rtm[selected], "latent_RTM": latent_rtm[selected]}
        print(f"Matched {sweep}: sensitivity={sensitivity:.1f}, FAR={far:.1f}; selected={check['n_selected']}, RTM={check['n_measured_RTM']}", flush=True)
    frame = pd.DataFrame([rows[key] for key in CONDITIONS])
    if frame.n_selected.tolist() != EXPECTED_SELECTED:
        raise AssertionError("Requested selection sizes differ from release.")
    if rows[(10,0)]["n_measurement_only"] or rows[(10,0)]["n_latent_only"]:
        raise AssertionError("Perfect control must have zero discordance.")
    audit = {field: np.concatenate([records[key][field] for key in CONDITIONS]) for field in records[CONDITIONS[0]]}
    summary_path, audit_path = output/"paired_summary.csv", output/"paired_selected_participants.npz"
    frame.to_csv(summary_path, index=False)
    np.savez_compressed(audit_path, **audit)
    mismatches = {k: {"release": v, "replay": software_record().get(k)} for k,v in release["software"].items() if software_record().get(k) != v}
    metadata = {"analysis": "Operational participant-level measurement-only RTM among detector-selected participants",
        "created_utc": datetime.now(timezone.utc).isoformat(), "source_commit": "59b47a5564342d4b2695307a341fbd0ff1d90091",
        "base_seed": release["base_seed"], "latent_seed": release["latent_seed"],
        "sensitivity_detection_seed": release["sensitivity_detection_seed"], "false_alarm_seed": release["false_alarm_seed"],
        "n_potential_participants": len(daily), "conditions": [list(x) for x in CONDITIONS],
        "condition_encoding": "Integer tenths: sensitivity, FAR per day; condition_id indexes this list; participant_id is original zero-based index.",
        "eligibility": "Each corrected baseline month >=3 and two-month mean >=4; no 25-day rule; no latent/perfect-detector eligibility filter.",
        "classification": "Strict B>reference and |T-reference|<B-reference, using exact integer comparisons; measured reference sensitivity*L/36; latent reference L/36.",
        "primary": "Measured RTM and not latent RTM divided by all selected in the same condition.",
        "secondary": "Measured RTM and not latent RTM divided by measured RTM in the same condition.",
        "uncertainty": "95% Wilson intervals from the unchanged canonical helper for measurement-only among selected and among measured RTM; conditional Monte Carlo/binomial uncertainty under this simulation, not real-world clinical uncertainty. No significance tests or between-condition confidence intervals.",
        "limitation": "An operational simulation-defined paired category covering measurement noise, FAR correction and flooring; not exclusive biological attribution or an additive causal decomposition. Gross discordance differs from net prevalence change and MPC. Within-condition pairing holds selection fixed; cross-condition differences combine selection and measurement. The independent 36-month reference is a finite-history estimate, not a known infinite-horizon mean.",
        "replayed_all_21_grid_rows_before_subsetting": True, "shared_perfect_control_checked": True,
        "latent_generation": latent_metadata, "software": software_record(), "environment_mismatches_from_release": mismatches,
        "match_scope": "Long-term total hash and selected/measured RTM counts at all 21 released grid rows match exactly; daily source hashes were not recorded by the release, so no claim of released daily-array byte identity is made.",
        "source_array_hashes": {k: array_record(a) for k,a in arrays.items()}, "detector_array_hashes": detector_hashes,
        "all_21_released_counts_match": all(x["exact_match"] for x in grid_checks),
        "independent_integer_classification_crosscheck": "All 100,000 participants in all 21 grid rows: exact agreement with canonical strict indicator for measured and latent classifications.",
        "source_file_hashes": {f: file_hash(ROOT/f) for f in ["realSim.py","reproduce_rtm3.py","rtm3_paired.py"]},
        "protected_released_file_hashes": released_hashes,
        "artifact_sha256": {p.name: file_hash(p) for p in [summary_path, audit_path, output/"replay_grid_validation.json"]},
        "participant_audit_arrays": {k: array_record(a) for k,a in audit.items()}}
    for path, sha in released_hashes.items():
        if file_hash(ROOT/path) != sha:
            raise AssertionError(f"Canonical artifact changed: {path}")
    (output/"paired_metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(frame.to_string(index=False), flush=True)
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/paired_type3")
    parser.add_argument("--latent-cache", type=Path, default=ROOT/"results/paired_type3/cache/latent_replay_cache.npz")
    args = parser.parse_args()
    run_paired_analysis(args.output_dir, args.latent_cache)


if __name__ == "__main__":
    main()
