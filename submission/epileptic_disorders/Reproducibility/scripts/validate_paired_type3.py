#!/usr/bin/env python3
"""Recompute all paired results from the compact participant audit, without replay.

Uses Python integer scalar arithmetic, not either production classification helper.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def validate(directory):
    directory = Path(directory)
    metadata = json.loads((directory/'paired_metadata.json').read_text())
    for name, expected in metadata['artifact_sha256'].items():
        assert hashlib.sha256((directory/name).read_bytes()).hexdigest() == expected, name
    audit = dict(np.load(directory/'paired_selected_participants.npz', allow_pickle=False))
    with (directory/'paired_summary.csv').open() as handle:
        summary = list(csv.DictReader(handle))
    independent=[]
    for condition, row in enumerate(summary):
        indices=np.flatnonzero(audit['condition_id'] == condition)
        counts=Counter(); diagnostics=Counter(); seen=set()
        for i in indices:
            ident=int(audit['participant_id'][i]); assert ident not in seen; seen.add(ident)
            assert 0 <= ident < metadata['n_potential_participants']
            true=list(map(int,audit['true_monthly'][i]))
            raw=list(map(int,audit['observed_monthly'][i]))
            measured=list(map(int,audit['measured_monthly'][i]))
            p=int(audit['sensitivity_tenths'][i]); q=10; far=int(audit['FAR_tenths'][i])
            assert [p,far] == metadata['conditions'][condition]
            assert measured == [max(x-3*far,0) for x in raw]
            assert min(measured[:2]) >= 3 and sum(measured[:2]) >= 8
            assert sum(measured[:2]) == audit['measured_baseline_total'][i]
            assert sum(measured[2:]) == audit['measured_test_total'][i]
            assert sum(true[:2]) == audit['true_baseline_total'][i]
            assert sum(true[2:]) == audit['true_test_total'][i]
            L=int(audit['long_term_total'][i])
            A,C=sum(measured[:2]),sum(measured[2:])
            m=(2*C < 3*A) and (9*q*A + 6*q*C > p*L)
            A,C=sum(true[:2]),sum(true[2:])
            l=(2*C < 3*A) and (9*A + 6*C > L)
            assert m == bool(audit['measured_RTM'][i])
            assert l == bool(audit['latent_RTM'][i])
            counts[(m,l)]+=1
            diagnostics['latent_eligible']+=min(true[:2])>=3 and sum(true[:2])>=8
            diagnostics['zero_true_baseline']+=sum(true[:2])==0
        n=len(indices); both=counts[(True,True)]; mo=counts[(True,False)]; lo=counts[(False,True)]; neither=counts[(False,False)]
        expected={'n_selected':n,'n_both_RTM':both,'n_measurement_only':mo,'n_latent_only':lo,'n_neither_RTM':neither,
                  'n_measured_RTM':both+mo,'n_latent_RTM':both+lo,
                  'n_selected_latent_eligible':diagnostics['latent_eligible'],
                  'n_selected_latent_ineligible':n-diagnostics['latent_eligible'],
                  'n_selected_zero_true_baseline':diagnostics['zero_true_baseline']}
        for key,value in expected.items(): assert int(row[key])==value,(condition,key)
        ratios={'measured_RTM_fraction':(both+mo)/n,'latent_RTM_fraction':(both+lo)/n,
                'measurement_only_fraction_selected':mo/n,'measurement_only_fraction_measured_RTM':mo/(both+mo),
                'latent_only_fraction_selected':lo/n,'net_measured_minus_latent_fraction':(mo-lo)/n}
        for key,value in ratios.items(): assert abs(float(row[key])-value) < 1e-15,(condition,key)
        for suffix,denominator in [('selected',n),('measured_RTM',both+mo)]:
            z=NormalDist().inv_cdf(.975); proportion=mo/denominator
            divider=1+z*z/denominator
            center=(proportion+z*z/(2*denominator))/divider
            half=z*math.sqrt(proportion*(1-proportion)/denominator+z*z/(4*denominator**2))/divider
            assert abs(float(row[f'measurement_only_{suffix}_ci_lo'])-(center-half)) < 1e-14
            assert abs(float(row[f'measurement_only_{suffix}_ci_hi'])-(center+half)) < 1e-14
        assert both+mo+lo+neither == n
        assert (both+mo)-(both+lo) == mo-lo
        independent.append({'condition_id':condition,**expected,'all_stored_fractions_match':True})
    for name,expected in metadata['protected_released_file_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    output={'status':'PASS','method':'Independent Python integer inequalities applied to every saved selected record; no production classification helper imported.',
            'records_checked':len(audit['participant_id']),'unique_conditions':5,
            'paired_cells_eligibility_correction_totals_reference_and_ratios_checked':True,
            'Wilson_intervals_independently_checked_with_NormalDist':True,
            'zero_true_baseline_and_latent_ineligible_retained':True,'canonical_files_unchanged':True,
            'conditions':independent,'participant_npz_sha256':metadata['artifact_sha256']['paired_selected_participants.npz']}
    (directory/'paired_validation.json').write_text(json.dumps(output,indent=2)+'\n')
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=ROOT/'results/paired_type3')
    print(json.dumps(validate(parser.parse_args().directory),indent=2))
