"""Scientific regression checks for the canonical RTM3 analysis."""
import hashlib
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import reproduce_rtm3 as rtm


class CoreSimulationTests(unittest.TestCase):
    def test_pinned_upstream_simulator_bytes(self):
        digest = hashlib.sha256((ROOT / "realSim.py").read_bytes()).hexdigest()
        self.assertEqual(digest, rtm.SIMULATOR_PROVENANCE["upstream_source_sha256"])

    def test_trial_and_reference_preserve_participant_parameters(self):
        np.random.seed(8181)
        statuses = set()
        for _ in range(200):
            trial, reference = rtm.generate_patient_histories(8.0)
            self.assertEqual(float(trial[1]), float(reference[1]))
            self.assertEqual(bool(trial[8]), bool(reference[8]))
            np.testing.assert_array_equal(trial[5], reference[5])
            np.testing.assert_array_equal(trial[6], reference[6])
            self.assertEqual(len(trial[0]), 150)
            self.assertEqual(len(reference[0]), 1080)
            statuses.add(bool(trial[8]))
        self.assertEqual(statuses, {False, True})

    def test_latent_generation_repeats_with_same_seed(self):
        first = rtm.generate_latent_cohort(30, 712, 0)
        second = rtm.generate_latent_cohort(30, 712, 0)
        for actual, expected in zip(first[:3], second[:3]):
            np.testing.assert_array_equal(actual, expected)
        self.assertEqual(first[0].dtype, np.int16)
        self.assertEqual(first[2].dtype, np.float64)
        self.assertTrue((first[0] >= 0).all())
        self.assertTrue((first[0] <= 144).all())

    def test_exact_rtm_distance_tie_is_false(self):
        monthly = np.array([[5, 5, 0, 1, 1], [5, 5, 0, 1, 2], [5, 5, 0, 0, 1]])
        reference = np.full(3, 102 / 36, dtype=np.float32)
        actual = rtm.observed_rtm_indicator(monthly, reference, 1.0)
        np.testing.assert_array_equal(actual, [False, True, False])

    def test_exact_endpoint_matches_integer_arithmetic_over_grid(self):
        rng = np.random.default_rng(729)
        monthly = rng.integers(0, 30, size=(1000, 5))
        long_total = rng.integers(0, 1200, size=1000)
        for sensitivity in rtm.SENSITIVITY_VALUES:
            base = 180 * monthly[:, :2].sum(axis=1)
            test = 120 * monthly[:, 2:].sum(axis=1)
            reference = int(round(10 * sensitivity)) * long_total
            expected = (base > reference) & (base - reference > np.abs(test - reference))
            actual = rtm.observed_rtm_indicator(monthly, long_total / 36, float(sensitivity), long_term_total=long_total)
            np.testing.assert_array_equal(actual, expected)

    def test_uncorrected_reference_includes_expected_alarms(self):
        monthly = np.array([[5, 5, 0, 1, 1], [5, 5, 0, 1, 2]])
        reference = np.full(2, 102 / 36)
        expected = rtm.observed_rtm_indicator(monthly, reference, 1.0)
        actual = rtm.observed_rtm_indicator(monthly + 30, reference, 1.0, reference_offset_monthly=30)
        np.testing.assert_array_equal(actual, expected)

    def test_fixed_cohort_reports_undefined_mpc_and_keeps_rtm_denominator(self):
        row = rtm.summarize_condition(
            monthly=np.array([[0, 0, 3, 3, 3], [4, 4, 3, 3, 3]]),
            sampled_monthly_rate=np.array([3., 3.]),
            effective_monthly_rate=np.array([3., 3.]),
            sensitivity=1., far=0., sweep="test", eligibility_variant="fixed",
            selection_mask=np.ones(2, dtype=bool), cohort_strategy=rtm.FIXED_COHORT,
        )
        self.assertEqual(row["n_eligible"], 2)
        self.assertEqual(row["n_valid_MPC"], 1)
        self.assertEqual(row["n_zero_baseline"], 1)
        self.assertEqual(row["n_undefined_MPC"], 1)
        self.assertEqual(row["MPC_valid_fraction_of_selected"], 0.5)
        self.assertEqual(row["MPC_median"], 25.)
        self.assertEqual(row["observed_RTM_fraction"], 0.5)

    def test_all_zero_fixed_cohort_has_no_defined_mpc(self):
        row = rtm.summarize_condition(
            monthly=np.zeros((2, 5)), sampled_monthly_rate=np.ones(2),
            effective_monthly_rate=np.ones(2), sensitivity=0.1, far=0.,
            sweep="test", eligibility_variant="fixed", selection_mask=np.ones(2, dtype=bool),
        )
        self.assertEqual(row["n_valid_MPC"], 0)
        self.assertTrue(np.isnan(row["MPC_median"]))
        self.assertEqual(row["observed_RTM_fraction"], 0.)

    def test_additive_thresholds_equal_corrected_eligibility(self):
        rng = np.random.default_rng(10)
        raw = rng.integers(0, 80, size=(5000, 5))
        for far in rtm.FAR_VALUES:
            expected = rtm.eligible_mask(np.maximum(raw - 30 * far, 0))
            actual = rtm.additive_threshold_mask(raw, float(far))
            np.testing.assert_array_equal(actual, expected)
        # Both individual months and their combined mean are tested at boundaries.
        boundary = np.array([[33, 35, 0, 0, 0], [32, 36, 0, 0, 0], [33, 34, 0, 0, 0]])
        np.testing.assert_array_equal(rtm.additive_threshold_mask(boundary, 1.), [True, False, False])

    def test_interval_accepts_25_but_rejects_26_zero_days(self):
        days = np.ones((2, 150), dtype=int)
        days[0, :25] = 0
        days[1, :26] = 0
        monthly = rtm.monthly_counts(days)
        np.testing.assert_array_equal(rtm.eligible_mask(monthly, days, True), [True, False])

    def test_small_median_sample_uses_unbounded_interval(self):
        median, lo, hi = rtm.exact_median_interval(np.array([5.]))
        self.assertEqual(median, 5.)
        self.assertEqual(lo, -np.inf)
        self.assertEqual(hi, np.inf)

    def test_sweeps_share_realizations_and_preserve_denominators(self):
        rng = np.random.default_rng(990)
        rates = rng.uniform(0.05, 2., size=500)
        daily = rng.poisson(rates[:, None], size=(500, 150)).astype(np.int16)
        effective = np.rint(rates * 30 * 36) / 36
        common = dict(latent_daily=daily, sampled_monthly_rate=rates * 30, effective_monthly_rate=effective)
        sens_primary, sens_interval, sens_cohort = rtm.run_sensitivity_sweep(**common, seed=991)
        far_primary, far_interval, correction, far_cohort, thresholds = rtm.run_far_sweep(**common, seed=992)
        primary = pd.DataFrame(sens_primary + far_primary)
        cohorts = pd.DataFrame(sens_cohort + far_cohort)
        self.assertEqual(len(cohorts), 42)
        self.assertEqual(len(thresholds), 22)
        summary = rtm.validate_results(
            primary, pd.DataFrame(sens_interval + far_interval), pd.DataFrame(correction),
            cohorts, pd.DataFrame(thresholds),
        )
        self.assertEqual(summary["status"], "PASS")
        reselected = cohorts[cohorts.cohort_strategy == rtm.RESELECTED_COHORT].reset_index(drop=True)
        pd.testing.assert_frame_equal(reselected, primary)
        raw_rows = pd.DataFrame(correction).query("correction_strategy == 'Uncorrected observed counts; 25-day rule omitted'")
        np.testing.assert_allclose(raw_rows.reference_offset_monthly, 30 * raw_rows.FAR)


if __name__ == "__main__":
    unittest.main()
