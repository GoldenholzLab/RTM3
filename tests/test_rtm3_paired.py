"""Boundary, reference-scale, selection and paired-cell tests for the addition."""
import unittest
import numpy as np
from rtm3_paired import classify_condition, paired_counts, strict_from_totals
from reproduce_rtm3 import observed_rtm_indicator


class PairedRTMTests(unittest.TestCase):
    def test_all_four_paired_cells(self):
        positive, negative = [4,4,2,2,2], [4,4,4,4,4]
        measured = np.array([positive,positive,negative,negative])
        latent = np.array([positive,negative,positive,negative])
        selected, m, l = classify_condition(latent, measured, np.full(4,72), 1.0)
        self.assertTrue(selected.all())
        row = paired_counts(m[selected],l[selected])
        self.assertEqual([row[k] for k in ['n_both_RTM','n_measurement_only','n_latent_only','n_neither_RTM']], [1,1,1,1])
        self.assertEqual(row['measurement_only_fraction_selected'],.25)
        self.assertEqual(row['measurement_only_fraction_measured_RTM'],.5)
        self.assertEqual(row['latent_only_fraction_selected'],.25)
        self.assertEqual(row['net_measured_minus_latent_fraction'],0)

    def test_exact_boundaries(self):
        # Reference=3. Equal upper/lower distances and baseline==reference fail.
        baseline = np.array([8,8,8,8,6,0])
        test = np.array([12,6,7,11,6,0])
        total = np.full(6,108)
        self.assertEqual(strict_from_totals(baseline,test,total).tolist(), [False,False,True,True,False,False])

    def test_reference_scaling(self):
        measured = np.array([[4,4,2,2,2]])
        total = np.array([180])
        selected,m,l = classify_condition(measured,measured,total,.5)
        self.assertTrue(selected[0])
        self.assertTrue(m[0])
        self.assertFalse(l[0])
        self.assertTrue(strict_from_totals(np.array([8]),np.array([6]),total,1,2)[0])
        self.assertFalse(strict_from_totals(np.array([8]),np.array([6]),total)[0])

    def test_exact_rational_reference_tie(self):
        # B=3, T=1, reference=2 from s=.1 and L=720: tied distance.
        monthly = np.array([[3,3,1,1,1]])
        total = np.array([720])
        self.assertFalse(observed_rtm_indicator(monthly,total/36,.1,long_term_total=total)[0])
        self.assertFalse(strict_from_totals(np.array([6]),np.array([3]),total,1,10)[0])

    def test_perfect_control_has_zero_discordance(self):
        monthly=np.array([[4,4,2,2,2],[4,4,4,4,4],[0,0,0,0,0]])
        selected,m,l=classify_condition(monthly,monthly,np.array([72,72,0]),1.0)
        self.assertTrue(np.array_equal(m,l))
        row=paired_counts(m[selected],l[selected])
        self.assertEqual(row['n_measurement_only'],0)
        self.assertEqual(row['n_latent_only'],0)

    def test_latent_zero_and_ineligible_are_retained(self):
        true=np.array([[0,0,0,0,0],[1,1,1,1,1]])
        measured=np.array([[4,4,1,1,1],[4,4,1,1,1]])
        selected,m,l=classify_condition(true,measured,np.array([36,36]),1.)
        self.assertEqual(selected.tolist(),[True,True])
        self.assertEqual(l.tolist(),[False,False])
        self.assertEqual(paired_counts(m[selected],l[selected])['n_selected'],2)

    def test_eligibility_inclusive_boundary(self):
        monthly=np.array([[3,5,1,1,1],[3,4,1,1,1],[2,6,1,1,1]])
        selected,_,_=classify_condition(monthly,monthly,np.full(3,36),1.)
        self.assertEqual(selected.tolist(),[True,False,False])

    def test_exhaustive_small_integer_crosscheck(self):
        # All 9,600 combinations; compare independent bounds with canonical rule.
        triples=np.array([(a,c,l) for a in range(12) for c in range(20) for l in range(40)])
        monthly=np.zeros((len(triples),5),dtype=np.int64)
        monthly[:,0]=triples[:,0]; monthly[:,2]=triples[:,1]
        for sensitivity, numerator in [(.1,1),(.5,5),(1.,10)]:
            expected=observed_rtm_indicator(monthly,triples[:,2]/36,sensitivity,long_term_total=triples[:,2])
            actual=strict_from_totals(triples[:,0],triples[:,1],triples[:,2],numerator,10)
            self.assertTrue(np.array_equal(actual,expected))


if __name__=='__main__':
    unittest.main()
