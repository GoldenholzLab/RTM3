#python code
import numpy as np
import pandas as pd
from realSim_turbo import simple_CHOCOLATES,downsample
import rtm_constants as CONST
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.auto import trange, tqdm
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from joblib import Parallel, delayed
import scipy.stats as stats
from itertools import groupby


def _require_baseline_eligibility(use_baselineTF):
    if not use_baselineTF:
        raise ValueError('Pre-baseline eligibility is retired for this paper; use_baselineTF must be True.')


# build one patient for trial as well as the recruitment period
def sim_one_potential_patient(recruitmentMonths, baselineMonths, testMonths, sensitivity, FAR, rng):
    # inputs
    #   recruitmentMonths - number of months for recruitment period
    #   baselineMonths - number of months for baseline period
    #   testMonths - number of months for test period
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    #   rng - random number generator to use 
    # outputs
    #   fullDiary - numpy array of daily detected seizure counts for the full patient period
    #   mSF - mean seizure frequency used in the simulation

    daysNeeded = (recruitmentMonths + baselineMonths + testMonths) * CONST.DAYS_PER_MONTH
    fullDiary, mSF = simple_CHOCOLATES(daysNeeded, defaultSeizureFreq=None , rng=rng)

    # apply sensitivity and FAR to full diary
    fullDeviceDiary = apply_sens_and_far_to_dairy(fullDiary, sensitivity, FAR, rng=rng)

    return fullDeviceDiary, mSF

# apply sensitivity and FAR to daily diary
def apply_sens_and_far_to_dairy(dailyDiary, sensitivity, FAR, rng):
    # inputs
    #   dailyDiary - numpy array of daily seizure counts
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    #   rng - random number generator to use 
    # outputs
    #   detectedDiary - numpy array of daily detected seizure counts
        # apply sensitivity probabilistically - each seizure has probability 'sensitivity' of being detected
    if sensitivity < 1.0:
        dailyDiary = np.array([rng.binomial(n=int(count), p=sensitivity) if count > 0 else 0 for count in dailyDiary])
    # apply FAR probabilistically    
    if FAR > 0.0:
        dailyDiary += rng.poisson(FAR,size=len(dailyDiary))

    return dailyDiary

def get_one_trial_patient(recruitmentMonths, baselineMonths, testMonths, sensitivity, FAR,
        eligibility_rate, min_monthly_rate, max_szfree_days, treatment_effect, use_baselineTF, 
    correct_FAR, rng, return_debug=False):
    # inputs
    #   recruitmentMonths - number of months for recruitment period
    #   baselineMonths - number of months for baseline period
    #   testMonths - number of months for test period
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    #   treatment_effect - reduction in monthly seizure frequency due to treatment (0 to 1)
    #   use_baselineTF - retained for compatibility; must be True so eligibility is tested during baseline
    #   correct_FAR - boolean, whether to correct for FAR in treatment effect
    #   rng - random number generator to use (this is optional, one will be made in CHOCOLATES if needed)
    # outputs   
    #   PC - percentage change in seizure frequency after treatment
    #   mSF - mean seizure frequency used in the simulation
    #   tries - number of attempts to find an eligible patient
    #   fails - number of times baseline period had zero seizures

    _require_baseline_eligibility(use_baselineTF)
    keepGoing = True
    tries = 0
    fails = 0
    debug_info = None
    PC = np.nan
    mSF = np.nan
    while keepGoing:
        fullDeviceDiary, mSF = sim_one_potential_patient(recruitmentMonths, baselineMonths, testMonths,
                                    sensitivity, FAR, rng)
        eligibleTF, baseline_start_idx = check_eligibility(fullDeviceDiary, recruitmentMonths, use_baselineTF,
                        correct_FAR, eligibility_rate, min_monthly_rate, max_szfree_days, FAR)
        if eligibleTF:
            PC, RTM_tf, BASEfailed, debug = treat_and_simplify_diary(
                fullDeviceDiary,
                baseline_start_idx,
                baselineMonths,
                testMonths,
                treatment_effect,
                correct_FAR,
                FAR,
                mSF,
                return_debug=return_debug,
            )
            if return_debug and debug is not None:
                debug['baseline_start_idx'] = int(baseline_start_idx)
                debug['RTM_tf'] = RTM_tf
                debug_info = debug
            keepGoing = BASEfailed  # if baseline failed, need to try again
            fails += BASEfailed # count how many times baseline failed
        else:
            keepGoing = True    # not eligible, try again
        tries += 1

    return PC, RTM_tf, mSF, tries, fails, debug_info

def treat_and_simplify_diary(fullDeviceDiary, baseline_start_idx,
    baselineMonths, testMonths, treatment_effect, correct_FAR, FAR, mSF, return_debug=False):
    # inputs
    #   fullDeviceDiary - numpy array of daily detected seizure counts for the full patient period
    #   baseline_start_idx - starting index of baseline period in fullDeviceDiary
    #   baselineMonths - number of months for baseline period
    #   testMonths - number of months for test period
    #   treatment_effect - reduction in monthly seizure frequency due to treatment (0 to 1)
    #   correct_FAR - boolean, whether to correct for FAR in treatment effect
    #   FAR - false alarm rate (false positives per day)
    #   mSF - monthly seizure frequency used in the simulation
    #   return_debug - boolean, whether to return debug information
    # outputs   
    #   PC - percentage change in seizure frequency after treatment
    #   BASEfailed - boolean, whether baseline period had zero seizures

    # extract baseline and test diaries
    baseline_end_idx = baseline_start_idx + baselineMonths * CONST.DAYS_PER_MONTH
    test_start_idx = baseline_end_idx
    test_end_idx = test_start_idx + testMonths * CONST.DAYS_PER_MONTH

    diaryBASE_daily = fullDeviceDiary[baseline_start_idx:baseline_end_idx]
    diaryTEST_daily = fullDeviceDiary[test_start_idx:test_end_idx]

    # downsample to monthly counts
    diaryBASE = downsample(diaryBASE_daily, CONST.DAYS_PER_MONTH)
    diaryTEST = downsample(diaryTEST_daily, CONST.DAYS_PER_MONTH)

    # note this tests RTM for this patient BEFORE applying drug.
    RTM_tf = check_RTM(mSF, diaryBASE, diaryTEST)

    # apply treatment effect
    diaryTEST = treat_one_patient_MONTHLY(diaryTEST, treatment_effect)

    # correct for known FAR if desired
    if correct_FAR and FAR > 0:
        mean_FAR_monthly = FAR * CONST.DAYS_PER_MONTH
        diaryBASE = diaryBASE - mean_FAR_monthly
        diaryBASE[diaryBASE < 0] = 0
        diaryTEST = diaryTEST - mean_FAR_monthly
        diaryTEST[diaryTEST < 0] = 0
    
    
    # compute percentage change in seizure frequency
    if np.sum(diaryBASE) > 0:
        PC = 100 * (np.mean(diaryBASE) - np.mean(diaryTEST)) / np.mean(diaryBASE)
        BASEfailed = False
    else:
        PC = 0
        BASEfailed = True

    debug = None
    if return_debug:
        debug = {
            'baseline_mean_monthly': float(np.mean(diaryBASE)),
            'test_mean_monthly': float(np.mean(diaryTEST)),
            'baseline_sum': float(np.sum(diaryBASE)),
            'test_sum': float(np.sum(diaryTEST)),
        }

    return PC, RTM_tf, BASEfailed, debug

def check_eligibility(fullDeviceDiary, recruitmentMonths, use_baselineTF,
            correct_FAR, eligibility_rate, min_monthly_rate, max_szfree_days, FAR):
    # inputs
    #   fullDeviceDiary - numpy array of daily detected seizure counts for the full patient period
    #   recruitmentMonths - number of months for recruitment period
    #   use_baselineTF - retained for compatibility; must be True so eligibility is tested during baseline
    #   correct_FAR - boolean, whether to correct for FAR in eligibility check
    #   eligibility_rate - fraction of mean monthly seizure frequency required for eligibility  
    #   min_monthly_rate - minimum monthly seizure frequency for eligibility
    #   max_szfree_days - maximum number of seizure-free days allowed in baseline for eligibility
    #   FAR - false alarm rate (false positives per day)
    #
    # outputs
    #   eligibleTF - boolean, whether the patient is eligible
    #   baseline_start_idx - starting index of baseline period in fullDeviceDiary
    _require_baseline_eligibility(use_baselineTF)
    search_start = 0
    search_end = recruitmentMonths * CONST.DAYS_PER_MONTH
    baseline_start_idx = np.nan
    eligibleTF = False

    for i in range(search_start, search_end - CONST.ELIGIBILITY_MONTHS * CONST.DAYS_PER_MONTH + 1):
        window = fullDeviceDiary[i:i + CONST.ELIGIBILITY_MONTHS * CONST.DAYS_PER_MONTH].copy()
        
        # Check average monthly rate
        monthly_counts = downsample(window, CONST.DAYS_PER_MONTH)
        if correct_FAR and FAR > 0:
            mean_FAR_monthly = FAR * CONST.DAYS_PER_MONTH
            monthly_counts = monthly_counts - mean_FAR_monthly
            monthly_counts[monthly_counts < 0] = 0

        avg_monthly_rate = np.mean(monthly_counts)
        
        if avg_monthly_rate < eligibility_rate:
            continue
        
        # Check minimum monthly rate
        if np.any(monthly_counts < min_monthly_rate):
            continue
        
        # Check for long seizure-free runs
        if correct_FAR and FAR > 0:
            window -= np.round(FAR).astype(int)
            window[window < 0] = 0.0
        zero_runs = (len(list(group)) for value, group in groupby(window) if value == 0)
        max_consecutive_zeros = max(zero_runs, default=0)
        if max_consecutive_zeros > max_szfree_days:
            continue
        
        # Criteria met
        eligibleTF = True
        if use_baselineTF:
            baseline_start_idx = i
        else:
            baseline_start_idx = i + CONST.ELIGIBILITY_MONTHS * CONST.DAYS_PER_MONTH
        break

    return eligibleTF, baseline_start_idx

# modify the test period seizure rate based on treatment effect
def treat_one_patient_MONTHLY(diaryTEST, treatment_effect):
    # inputs
    #   diaryTEST - numpy array of monthly seizure counts
    #   treatment_effect - reduction in monthly seizure frequency due to treatment
    # outputs
    #   treated_diaryTEST - numpy array of monthly seizure counts after treatment

    # if there is no effect, do nothing
    if treatment_effect <= 0:
        return diaryTEST
    
    treated_diaryTEST = diaryTEST * (1 - treatment_effect)
    treated_diaryTEST[treated_diaryTEST < 0] = 0  # ensure no negative seizure counts
    return treated_diaryTEST

def build_one_trial(numPatients, recruitmentMonths, baselineMonths, testMonths, sensitivity, FAR,
        eligibility_rate, min_monthly_rate, max_szfree_days, treatment_effect, use_baselineTF, 
        correct_FAR):
    # inputs
    #   numPatients - number of patients to simulate in the trial
    #   recruitmentMonths - number of months for recruitment period
    #   baselineMonths - number of months for baseline period
    #   testMonths - number of months for test period
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    #   treatment_effect - reduction in monthly seizure frequency due to treatment (0 to 1)
    #   use_baselineTF - retained for compatibility; must be True so eligibility is tested during baseline
    #   correct_FAR - boolean, whether to correct for FAR in treatment effect
    # outputs   
    #   PCs - list of percentage changes in seizure frequency after treatment for each patient
    #   mSFs - list of mean seizure frequencies used in the simulation for each patient
    #   total_tries - total number of attempts to find eligible patients
    #   total_fails - total number of times baseline period had zero seizures

    _require_baseline_eligibility(use_baselineTF)
    rng = np.random.default_rng()
    drug_PCs, placebo_PCs, mSFs = [], [], []
    total_tries, total_fails, rtm_count = 0, 0, 0
    halfPatients = numPatients // 2

    for ptNum in range(numPatients):
        if ptNum >= halfPatients:
            this_treatment_effect = 0    # no treatment effect for control group
        else:
            this_treatment_effect = treatment_effect  # treatment effect for treatment group

        PC, RTM_tf, mSF, tries, fails, _debug = get_one_trial_patient(
            recruitmentMonths,
            baselineMonths,
            testMonths,
            sensitivity,
            FAR,
            eligibility_rate,
            min_monthly_rate,
            max_szfree_days,
            this_treatment_effect,
            use_baselineTF,
            correct_FAR,
            rng,
        )
        if this_treatment_effect > 0:    #treatment group
            drug_PCs.append(PC)
        else:                           #control group
            placebo_PCs.append(PC)
        mSFs.append(mSF)
        total_tries += tries
        total_fails += fails
        rtm_count += RTM_tf

    drug_PCs = np.array(drug_PCs)
    placebo_PCs =  np.array(placebo_PCs)

    results = {
        'drug_MPC': np.median(drug_PCs),
        'placebo_MPC': np.median(placebo_PCs),
        'drug_RR50': 100*np.mean(drug_PCs >= 50),
        'placebo_RR50': 100*np.mean(placebo_PCs >= 50),
        'success_RR50': calculate_fisher_exact_p_value(placebo_PCs, drug_PCs)<0.05,
        'success_MPC': calculate_MPC_p_value(placebo_PCs, drug_PCs)<0.05,
        'eligibility_rate': 100*numPatients / (total_tries + numPatients),
        'failed_baseline_rate': 100*total_fails / (total_fails + numPatients),
        'median_mSF': np.median(mSFs),
        'fraction_RTM': rtm_count / numPatients
    }
    return results

def check_RTM(mSF, monthlyBASELINE, monthlyTEST):
    mean_BASELINE = np.mean(monthlyBASELINE)
    if mean_BASELINE > mSF:
        if (mean_BASELINE - mSF) > np.abs(np.mean(monthlyTEST) - mSF):
            return True
    return False

def calculate_MPC_p_value(placebo_arm_percent_changes,drug_arm_percent_changes):

    # Mann_Whitney_U test
    [_, MPC_p_value] = stats.ranksums(placebo_arm_percent_changes, drug_arm_percent_changes)

    return MPC_p_value

def calculate_fisher_exact_p_value(placebo_arm_percent_changes,
                                drug_arm_percent_changes):

    num_placebo_arm_responders     = np.sum(placebo_arm_percent_changes >= 50)
    num_drug_arm_responders        = np.sum(drug_arm_percent_changes    >= 50)
    num_placebo_arm_non_responders = len(placebo_arm_percent_changes) - num_placebo_arm_responders
    num_drug_arm_non_responders    = len(drug_arm_percent_changes)    - num_drug_arm_responders

    table = np.array([[num_placebo_arm_responders, num_placebo_arm_non_responders], [num_drug_arm_responders, num_drug_arm_non_responders]])

    [_, RR50_p_value] = stats.fisher_exact(table)

    return RR50_p_value

def run_a_set_of_trials(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
        sensitivity, FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
        treatment_effect, use_baselineTF, correct_FAR):
    # inputs
    #   numTrials - number of trials to simulate
    #   numPatients - number of patients to simulate in each trial
    #   recruitmentMonths - number of months for recruitment period
    #   baselineMonths - number of months for baseline period
    #   testMonths - number of months for test period
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    #   treatment_effect - reduction in monthly seizure frequency due to treatment (0 to 1)
    #   use_baselineTF - retained for compatibility; must be True so eligibility is tested during baseline
    #   correct_FAR - boolean, whether to correct for FAR in treatment effect
    # outputs:   
    #   all_results - list of dictionaries containing analysis results for each trial
    _require_baseline_eligibility(use_baselineTF)

    all_results = Parallel(n_jobs=25)(
        delayed(build_one_trial)(
            numPatients, recruitmentMonths, baselineMonths, testMonths,
            sensitivity, FAR, eligibility_rate, min_monthly_rate,
            max_szfree_days, treatment_effect, use_baselineTF,
            correct_FAR
        )
        for _ in trange(numTrials, desc='Simulating Trials')
    )

    summary = {
        'numPatients': numPatients,
        'numTrials': numTrials,
        'recruitmentMonths': recruitmentMonths,
        'baselineMonths': baselineMonths,
        'testMonths': testMonths,
        'sensitivity': sensitivity,
        'FAR': FAR,
        'eligibility_rate': eligibility_rate,
        'min_monthly_rate': min_monthly_rate,
        'max_szfree_days': max_szfree_days,
        'treatment_effect': treatment_effect,
        'use_baselineTF': use_baselineTF,
        'correct_FAR': correct_FAR,
        'mean_drugMPC': np.mean([res['drug_MPC'] for res in all_results]),
        'mean_placeboMPC': np.mean([res['placebo_MPC'] for res in all_results]),
        'mean_diffMPC': np.mean([res['drug_MPC'] - res['placebo_MPC'] for res in all_results]),
        'mean_drugRR50': np.mean([res['drug_RR50'] for res in all_results]),
        'mean_placeboRR50': np.mean([res['placebo_RR50'] for res in all_results]),
        'mean_diffRR50': np.mean([res['drug_RR50'] - res['placebo_RR50'] for res in all_results]),
        'successRR50': 100*np.mean([res['success_RR50'] for res in all_results]),
        'successMPC': 100*np.mean([res['success_MPC'] for res in all_results]),
        'mean_mSF': np.mean([res['median_mSF'] for res in all_results]),
        'mean_eligibility_rate': np.mean([res['eligibility_rate'] for res in all_results]),
        'mean_failed_baseline_rate': np.mean([res['failed_baseline_rate'] for res in all_results]),
        'mean_fraction_RTM': np.mean([res['fraction_RTM'] for res in all_results]),
    }
    return pd.DataFrame(all_results), pd.DataFrame([summary])

# baseline-only eligibility summary
def test1(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
        sensitivity, FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
        treatment_effect, correct_FAR, fn):

    for i, use_baselineTF in enumerate([True]):
        allresults,summary = run_a_set_of_trials(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
            sensitivity, FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
            treatment_effect, use_baselineTF, correct_FAR)
        print(summary)
        if i==0:
            summary.to_csv(fn, index=False, mode='w', header=True)
        else:  
            summary.to_csv(fn, index=False, mode='a', header=False)
    return

# effect of sensitivity
def test2(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
        FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
        treatment_effect, correct_FAR, use_baselineTF, fn):

    for i, sensitivity in enumerate([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]):
        allresults,summary = run_a_set_of_trials(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
            sensitivity, FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
            treatment_effect, use_baselineTF, correct_FAR)
        print(summary)
        if i==0:
            summary.to_csv(fn, index=False, mode='w', header=True)
        else:  
            summary.to_csv(fn, index=False, mode='a', header=False)
    return

# check FAR
def test3(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
        sensitivity, eligibility_rate, min_monthly_rate, max_szfree_days,
        treatment_effect, use_baselineTF, fn):

    for i, FAR in enumerate([0,1/30,1/7,1]):
        for k, correct_FAR in enumerate([False, True]):
            allresults,summary = run_a_set_of_trials(numTrials, numPatients, recruitmentMonths, baselineMonths, testMonths,
                sensitivity, FAR, eligibility_rate, min_monthly_rate, max_szfree_days,
                treatment_effect, use_baselineTF, correct_FAR)
            print(summary)
            if i==0 and k==0:
                summary.to_csv(fn, index=False, mode='w', header=True)
            else:  
                summary.to_csv(fn, index=False, mode='a', header=False)
    return

def diagnostic_baselineTF_patient_level(
    num_patients_per_arm,
    recruitmentMonths,
    baselineMonths,
    testMonths,
    sensitivity,
    FAR,
    eligibility_rate,
    min_monthly_rate,
    max_szfree_days,
    treatment_effect,
    correct_FAR,
    seed=0,
):
    rng = np.random.default_rng(seed)
    rows = []

    for use_baselineTF in [True]:
        for arm, this_treatment_effect in [('drug', treatment_effect), ('placebo', 0.0)]:
            for _ in trange(
                num_patients_per_arm,
                desc=f'Patients baselineTF={use_baselineTF} arm={arm}',
                leave=False,
            ):
                PC, RTM_tf, mSF, tries, fails, debug = get_one_trial_patient(
                    recruitmentMonths,
                    baselineMonths,
                    testMonths,
                    sensitivity,
                    FAR,
                    eligibility_rate,
                    min_monthly_rate,
                    max_szfree_days,
                    this_treatment_effect,
                    use_baselineTF,
                    correct_FAR,
                    rng,
                    return_debug=True,
                )
                rows.append(
                    {
                        'use_baselineTF': use_baselineTF,
                        'arm': arm,
                        'PC': PC,
                        'baseline_mean_monthly': debug['baseline_mean_monthly'] if debug else np.nan,
                        'test_mean_monthly': debug['test_mean_monthly'] if debug else np.nan,
                        'baseline_sum': debug['baseline_sum'] if debug else np.nan,
                        'test_sum': debug['test_sum'] if debug else np.nan,
                        'baseline_start_idx': debug['baseline_start_idx'] if debug else np.nan,
                        'RTM_tf': RTM_tf,
                        'mSF': mSF,
                        'tries': tries,
                        'fails': fails,
                    }
                )

    df = pd.DataFrame(rows)
    summary = (
        df.groupby(['use_baselineTF', 'arm'])
        .agg(
            n=('PC', 'size'),
            PC_mean=('PC', 'mean'),
            PC_median=('PC', 'median'),
            PC_sd=('PC', 'std'),
            baseline_mean_mean=('baseline_mean_monthly', 'mean'),
            test_mean_mean=('test_mean_monthly', 'mean'),
            tries_mean=('tries', 'mean'),
            fails_rate=('fails', 'mean'),
        )
        .reset_index()
    )

    placebo = df[df['arm'] == 'placebo']
    placebo_q = (
        placebo.groupby('use_baselineTF')['PC']
        .quantile([0.05, 0.25, 0.5, 0.75, 0.95])
        .unstack()
        .reset_index()
        .rename(columns={0.05: 'q05', 0.25: 'q25', 0.5: 'q50', 0.75: 'q75', 0.95: 'q95'})
    )

    return df, summary, placebo_q
