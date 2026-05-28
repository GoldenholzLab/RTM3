#python code
import numpy as np
import pandas as pd
from realSim import get_mSF,simulator_base,downsample
import rtm_constants as CONST
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.auto import trange, tqdm
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from joblib import Parallel, delayed
import scipy.stats as stats
from plotting_tools import analyze_placebo_response


def _require_baseline_eligibility(use_baseline):
    if not use_baseline:
        raise ValueError('Pre-baseline eligibility is retired for this paper; use_baseline must be True.')


# modify the test period seizure rate based on treatment effect
def treat_one_patient(diaryTEST, treatment_effect):
    # inputs
    #   diaryTEST - numpy array of monthly seizure counts
    #   treatment_effect - reduction in monthly seizure frequency due to treatment
    # outputs
    #   treated_diaryTEST - numpy array of monthly seizure counts after treatment

    treated_diaryTEST = diaryTEST * (1 - treatment_effect)
    if np.isscalar(treated_diaryTEST) or np.ndim(treated_diaryTEST) == 0:
        return max(float(treated_diaryTEST), 0.0)
    treated_diaryTEST[treated_diaryTEST < 0] = 0  # ensure no negative seizure counts
    return treated_diaryTEST

# choose a reasonable monthly seizure frequency (wrapper)
# uses CHCOLATES to select from a log-noraml distribution
# such that realistic seizure frequencies are chosen
def choose_seizure_frequency():
    return get_mSF(requested_msf=-1)


def simplified_patient_generator(mSF=None, scale=CONST.MONTHLY_STD):
    if mSF is None:
        mSF = choose_seizure_frequency()

    diaryPRE = np.random.normal(loc=mSF, scale=scale, size=1)
    diaryBASELINE = np.random.normal(loc=mSF, scale=scale, size=1)
    diaryTEST = np.random.normal(loc=mSF, scale=scale, size=1)

    diaryPRE[diaryPRE < 0] = 0
    diaryBASELINE[diaryBASELINE < 0] = CONST.EPSILON
    diaryTEST[diaryTEST < 0] = 0

    return diaryPRE, diaryBASELINE, diaryTEST, mSF


def keepEligible(cohort_data, eligibility_min, use_baseline=True):
    _require_baseline_eligibility(use_baseline)
    eligible_data = []
    for patient in cohort_data:
        if patient['diaryBASELINE'] >= eligibility_min:
            eligible_data.append(patient)
    return eligible_data


def generate_simplified_cohort(num_patients=CONST.COHORT_SIZE, STD=CONST.MONTHLY_STD, eligibility_min=0):
    cohort_data = []
    per_arm = num_patients // 2

    for _ in range(per_arm):
        diaryPRE, diaryBASELINE, diaryTEST, mSF = simplified_patient_generator(mSF=None, scale=STD)
        cohort_data.append({
            'diaryPRE': diaryPRE[0],
            'diaryBASELINE': diaryBASELINE[0],
            'diaryTEST': diaryTEST[0],
            'mSF': mSF[0],
            'PC': 100*(1 - diaryTEST[0] / diaryBASELINE[0])
        })

    return pd.DataFrame(keepEligible(cohort_data, eligibility_min, use_baseline=True))


def generated_one_complex_patient(use_baseline=True, sensitivity=1.0, FAR=0.0):
    _require_baseline_eligibility(use_baseline)
    preM = CONST.ELIGIBILITY_MONTHS - CONST.TRIAL_BASELINE_MONTHS
    daysNeeded = (preM + CONST.TRIAL_BASELINE_MONTHS + CONST.TRIAL_TEST_MONTHS) * CONST.DAYS_PER_MONTH
    mSF = choose_seizure_frequency()

    dailyDiary = simulator_base(sampRATE=1, number_of_days=daysNeeded, defaultSeizureFreq=mSF)
    monthlyDiary = downsample(x=dailyDiary, byHowmuch=CONST.DAYS_PER_MONTH)

    monthlyDiary = (sensitivity * monthlyDiary).astype(int)
    monthlyDiary = monthlyDiary + np.random.poisson(FAR*CONST.DAYS_PER_MONTH, size=len(monthlyDiary))
    monthlyDiary = monthlyDiary - int(FAR*CONST.DAYS_PER_MONTH)
    monthlyDiary[monthlyDiary < 0] = 0

    diaryPre = np.mean(monthlyDiary[:preM]) if preM > 0 else np.nan
    diaryBaseline = np.mean(monthlyDiary[preM:(preM + CONST.TRIAL_BASELINE_MONTHS)])
    if diaryBaseline < CONST.EPSILON:
        diaryBaseline = CONST.EPSILON
    diaryTest = np.mean(monthlyDiary[(preM + CONST.TRIAL_BASELINE_MONTHS):])

    # compute an effective mSF
    dailyDiary = simulator_base(sampRATE=1, number_of_days=CONST.DAYS_PER_MONTH*36, defaultSeizureFreq=mSF)
    monthlyDiary36 = downsample(x=dailyDiary, byHowmuch=CONST.DAYS_PER_MONTH)
    effective_mSF = np.mean(monthlyDiary36)

    return diaryPre, diaryBaseline, diaryTest, effective_mSF


def generate_complex_cohort(num_patients=CONST.COHORT_SIZE, eligibility_min=0, use_baseline=True,
                            sensitivity=1.0, FAR=0.0):
    _require_baseline_eligibility(use_baseline)
    cohort_data = []
    per_arm = num_patients // 2
    for _ in range(per_arm):
        diaryPre, diaryBASELINE, diaryTEST, mSF = generated_one_complex_patient(
            use_baseline=use_baseline, sensitivity=sensitivity, FAR=FAR)
        cohort_data.append({
            'diaryPRE': diaryPre,
            'diaryBASELINE': diaryBASELINE,
            'diaryTEST': diaryTEST,
            'mSF': mSF, # this used to be mSF[0] but mSF is now a scalar, not an array
            'PC': 100*(1 - diaryTEST / diaryBASELINE)
        })
    return pd.DataFrame(keepEligible(cohort_data, eligibility_min, use_baseline=use_baseline))


def generated_one_complex_patient_v2(use_baseline, 
            sensitivity, FAR, eligibility_min, eligibility_monthly_min,
            eligibility_longest_szfree,correct_FAR):
    _require_baseline_eligibility(use_baseline)
    if use_baseline:
        pre_months = 0
        pre_days = 0
        eligibility_days = CONST.TRIAL_BASELINE_MONTHS * CONST.DAYS_PER_MONTH
    else:
        pre_months = CONST.ELIGIBILITY_MONTHS
        pre_days = pre_months * CONST.DAYS_PER_MONTH
        eligibility_days = pre_days

    daysNeeded = pre_days + (CONST.TRIAL_BASELINE_MONTHS + CONST.TRIAL_TEST_MONTHS) * CONST.DAYS_PER_MONTH

    mSF = choose_seizure_frequency()
    dailyDiary= simulator_base(sampRATE=1, number_of_days=daysNeeded, defaultSeizureFreq=mSF)
    
    # apply sensitivity probabilistically - each seizure has probability 'sensitivity' of being detected
    if sensitivity < 1.0:
        dailyDiary = np.array([np.random.binomial(n=int(count), p=sensitivity) if count > 0 else 0 for count in dailyDiary])
    # apply FAR probabilistically    
    if FAR > 0.0:
        dailyDiary += np.random.poisson(FAR,size=len(dailyDiary))
        
    # downsample to monthly    
    monthlyDiary = downsample(x=dailyDiary,byHowmuch=CONST.DAYS_PER_MONTH)
    if correct_FAR:
        # adjust for FAR
        remover = np.round(FAR * CONST.DAYS_PER_MONTH).astype(int)
        monthlyDiary = np.array([max(count - remover, 0.0) for count in monthlyDiary])

    eligibleTF = check_eligibility(dailyDiary[:eligibility_days], 
                monthlyDiary[:(eligibility_days//CONST.DAYS_PER_MONTH)],
                eligibility_min=eligibility_min,eligibility_monthly_min=eligibility_monthly_min,
                eligibility_longest_szfree=eligibility_longest_szfree,correct_FAR=correct_FAR, FAR=FAR)

    baseline_start = pre_months
    baseline_end = baseline_start + CONST.TRIAL_BASELINE_MONTHS
    diaryBaseline = np.mean(monthlyDiary[baseline_start:baseline_end])
    if diaryBaseline < CONST.EPSILON:
        diaryBaseline = CONST.EPSILON  # avoid divide by zero

    test_start = baseline_end
    test_end = test_start + CONST.TRIAL_TEST_MONTHS
    diaryTest = np.mean(monthlyDiary[test_start:test_end])    

    # compute an effective mSF
    dailyDiary = simulator_base(sampRATE=1, number_of_days=CONST.DAYS_PER_MONTH*36, defaultSeizureFreq=mSF)
    monthlyDiary36 = downsample(x=dailyDiary, byHowmuch=CONST.DAYS_PER_MONTH)
    effective_mSF = np.mean(monthlyDiary36)

    return eligibleTF, diaryBaseline, diaryTest, effective_mSF


def check_eligibility(dailyDiary, monthlyDiary, eligibility_min, eligibility_monthly_min, 
            eligibility_longest_szfree, correct_FAR, FAR):

    
    ### test monthly rate first
    baseline_avg = np.mean(monthlyDiary)
    
    if baseline_avg < eligibility_min:
        return False
    
    ### test each monthly minimum
    if np.any(monthlyDiary < eligibility_monthly_min):
        return False
    
    ### test longest seizure-free period
    if correct_FAR:
        dailyDiary -= np.round(FAR).astype(int)
        dailyDiary[dailyDiary < 0] = 0.0

    max_szfree = 0
    current_szfree = 0
    for day_count in dailyDiary:
        if day_count == 0:
            current_szfree += 1
            if current_szfree > max_szfree:
                max_szfree = current_szfree
        else:
            current_szfree = 0
    if max_szfree > eligibility_longest_szfree:
        return False
    
    # if all tests passed
    return True


def generate_complex_cohort_v2(num_patients,eligibility_min,use_baseline,
                            sensitivity,FAR,
                            eligibility_monthly_min,
                            eligibility_longest_szfree,
                            correct_FAR):
    _require_baseline_eligibility(use_baseline)

    cohort_data = []
    for _ in range(num_patients):
        eligibleTF, diaryBASELINE, diaryTEST, mSF = generated_one_complex_patient_v2(
            use_baseline=use_baseline, sensitivity=sensitivity, FAR=FAR,
            eligibility_min=eligibility_min,
            eligibility_monthly_min=eligibility_monthly_min,
            eligibility_longest_szfree=eligibility_longest_szfree,
            correct_FAR=correct_FAR)
        cohort_data.append({
                    'eligibleTF': eligibleTF,
                    'RTM_tf': RTM_tester(diaryBASELINE, diaryTEST, mSF, sensitivity),
                    'diaryBASELINE': diaryBASELINE,
                    'diaryTEST': diaryTEST,
                    'mSF': np.asarray(mSF).reshape(-1)[0],
                    'PC': 100*(1 - diaryTEST / diaryBASELINE)
                })
    return pd.DataFrame(cohort_data)    

def RTM_tester(diaryBASELINE, diaryTEST, mSF_real, sensitivity):
    mSF = mSF_real * sensitivity
    if diaryBASELINE > mSF:
        if (diaryBASELINE - mSF) > np.abs(diaryTEST - mSF):
            return True
    return False

def process_sensitivity(sensitivity,FAR,num_patients,correct_FAR,use_baseline,sz_min=4,sz_monthly_min=3,sz_longest_free=25):
    _require_baseline_eligibility(use_baseline)
    cohort_df = generate_complex_cohort_v2(num_patients=num_patients,
                        eligibility_min=sz_min, use_baseline=use_baseline,
                        sensitivity=sensitivity, FAR=FAR,
                        eligibility_monthly_min=sz_monthly_min,
                        eligibility_longest_szfree=sz_longest_free,
                        correct_FAR=correct_FAR)
    fraction_eligible = np.sum(cohort_df['eligibleTF']) / num_patients
    cohort2 = cohort_df[cohort_df['eligibleTF'] == True]
    if len(cohort2) == 0:
        frac_RTM = np.nan
        MPC = np.nan
    else:
        frac_RTM = np.sum(cohort2['RTM_tf']) / len(cohort2)
        MPC = np.nanmedian(cohort2['PC'])
    return {
        'sensitivity': sensitivity,
        'FAR': FAR,
        'use_baseline': use_baseline,
        'correct_FAR': correct_FAR,
        'fraction_eligible': fraction_eligible,
        'frac_RTM': frac_RTM,
        'MPC': MPC,
        'ekligibility_min': sz_min,
        'eligibility_monthly_min': sz_monthly_min,
        'eligibility_longest_szfree': sz_longest_free
    }

def _median_with_ci(values, confidence=0.95):
    values = np.asarray(values, dtype=float)
    values = np.sort(values[~np.isnan(values)])
    n_values = len(values)
    if n_values == 0:
        return np.nan, np.nan, np.nan
    if n_values == 1:
        median = float(values[0])
        return median, median, median

    alpha = 1.0 - confidence
    lower_index = int(stats.binom.ppf(alpha / 2.0, n_values, 0.5)) - 1
    upper_index = int(stats.binom.isf(alpha / 2.0, n_values, 0.5))
    lower_index = max(0, lower_index)
    upper_index = min(n_values - 1, upper_index)
    return float(np.nanmedian(values)), float(values[lower_index]), float(values[upper_index])

def process_sensitivity_mpc_ci(sensitivity,FAR,num_patients,correct_FAR,use_baseline,
                               sz_min=4,sz_monthly_min=3,sz_longest_free=25,
                               confidence=0.95):
    _require_baseline_eligibility(use_baseline)
    cohort_df = generate_complex_cohort_v2(num_patients=num_patients,
                        eligibility_min=sz_min, use_baseline=use_baseline,
                        sensitivity=sensitivity, FAR=FAR,
                        eligibility_monthly_min=sz_monthly_min,
                        eligibility_longest_szfree=sz_longest_free,
                        correct_FAR=correct_FAR)
    fraction_eligible = np.sum(cohort_df['eligibleTF']) / num_patients
    cohort2 = cohort_df[cohort_df['eligibleTF'] == True]
    if len(cohort2) == 0:
        frac_RTM = np.nan
        mpc_median, mpc_ci_lo, mpc_ci_hi = np.nan, np.nan, np.nan
    else:
        frac_RTM = np.sum(cohort2['RTM_tf']) / len(cohort2)
        mpc_median, mpc_ci_lo, mpc_ci_hi = _median_with_ci(cohort2['PC'], confidence=confidence)
    return {
        'sensitivity': sensitivity,
        'FAR': FAR,
        'use_baseline': use_baseline,
        'correct_FAR': correct_FAR,
        'fraction_eligible': fraction_eligible,
        'n_eligible': len(cohort2),
        'frac_RTM': frac_RTM,
        'MPC': mpc_median,
        'MPC_median': mpc_median,
        'MPC_ci_lo': mpc_ci_lo,
        'MPC_ci_hi': mpc_ci_hi,
        'MPC_ci_confidence': confidence,
        'ekligibility_min': sz_min,
        'eligibility_monthly_min': sz_monthly_min,
        'eligibility_longest_szfree': sz_longest_free
    }

def run_mpc_ci_sensitivity_and_far(use_baseline=True,num_patients=50000,fn='rtm_test123_mpc_ci_results.csv',
                                   confidence=0.95,sensitivity_values=None,far_values=None):
    _require_baseline_eligibility(use_baseline)
    if sensitivity_values is None:
        sensitivity_values = np.linspace(0.1, 1.0, 10)
    if far_values is None:
        far_values = np.linspace(0.0, 1.0, 11)

    rows = []
    for correct_FAR in (True, False):
        for sensitivity in sensitivity_values:
            rows.append(process_sensitivity_mpc_ci(
                sensitivity, 0.0, num_patients, correct_FAR, use_baseline,
                sz_min=4, sz_monthly_min=3, sz_longest_free=25,
                confidence=confidence))
    for correct_FAR in (True, False):
        for FAR in far_values:
            if FAR == 0.0:
                continue
            rows.append(process_sensitivity_mpc_ci(
                1.0, FAR, num_patients, correct_FAR, use_baseline,
                sz_min=4, sz_monthly_min=3, sz_longest_free=25,
                confidence=confidence))

    results = pd.DataFrame(rows)
    results.to_csv(fn, index=False)
    return results

def test1_sensitivity_par(use_baseline=True,num_patients=1000000):
    # note the eligibility rules come from Krauss et al 2022 cenobamate trial
    _require_baseline_eligibility(use_baseline)
    results = []
    correct_FAR = True
    results = Parallel(n_jobs=-1)(delayed(process_sensitivity)
                (s, 0.0,num_patients,correct_FAR,use_baseline,sz_min=4,sz_monthly_min=3,sz_longest_free=25) for s in np.linspace(0.1, 1.0, 10))
    return pd.DataFrame(results)

def test2_FAR_par(use_baseline=True, correct_FAR=True,num_patients=1000000):
    # note the eligibility rules come from Krauss et al 2022 cenobamate trial
    _require_baseline_eligibility(use_baseline)

    results = []
    results = Parallel(n_jobs=-1)(delayed(process_sensitivity)
                (1.0, f,num_patients,correct_FAR,use_baseline,sz_min=4,sz_monthly_min=3,sz_longest_free=25) for f in np.linspace(0.0, 1.0, 11))
    return pd.DataFrame(results)


def process_efficacy(sensitivity,FAR,num_patients,correct_FAR,use_baseline,
                sz_min,sz_monthly_min,sz_longest_free, drug_effect, Nlist):
    _require_baseline_eligibility(use_baseline)

    ## generate until we have enough eligible patients
    pcb_PC = []
    drug_PC = []
    mSF_list = []
    rtm_list = []
    cohortCount = 0
    rejectCount = 0
    drugCount = 0
    
    while cohortCount < num_patients:
        eligibleTF, diaryBASELINE, diaryTEST, mSF = generated_one_complex_patient_v2(
            use_baseline=use_baseline, sensitivity=sensitivity, FAR=FAR,
            eligibility_min=sz_min,
            eligibility_monthly_min=sz_monthly_min,
            eligibility_longest_szfree=sz_longest_free,
            correct_FAR=correct_FAR)
        if eligibleTF:    
            cohortCount += 1
            mSF_list.append(mSF)
            rtm_list.append(RTM_tester(diaryBASELINE, diaryTEST, mSF, sensitivity))
            if drugCount< (num_patients // 2):
                # treat this patient
                diaryTEST = treat_one_patient(diaryTEST, treatment_effect=drug_effect)
                drugCount += 1
                drug_PC.append(100*(1 - diaryTEST / diaryBASELINE))
            else:
                pcb_PC.append(100*(1 - diaryTEST / diaryBASELINE))
        else:
            rejectCount += 1
    
    fraction_eligible = (num_patients) / (num_patients + rejectCount)       # not using but interesting for $ analysis

    ### calcualte success of trial
    success= np.zeros(len(Nlist))
    successR= np.zeros(len(Nlist))
    for i,N in enumerate(Nlist):
        pcbN_PC = np.random.choice(pcb_PC, size=N, replace=True)
        drugN_PC = np.random.choice(drug_PC, size=N, replace=True)
        success[i]  = (calculate_MPC_p_value(pcbN_PC, drugN_PC) < 0.05)
        successR[i] = (calculate_fisher_exact_p_value(pcbN_PC, drugN_PC) < 0.05)
    ##
    MPC_drug_maxN = np.nanmedian(drug_PC) 
    MPC_pcb_maxN = np.nanmedian(pcb_PC)
    result = {
        'sensitivity': sensitivity,
        'FAR': FAR,
        'use_baseline': use_baseline,
        'correct_FAR': correct_FAR,
        'sz_min': sz_min,
        'sz_monthly_min': sz_monthly_min,
        'sz_longest_free': sz_longest_free,
        'fraction_eligible': fraction_eligible,
        'num_patients': num_patients,
        'MPC_drug_maxN': MPC_drug_maxN,
        'MPC_pcb_maxN': MPC_pcb_maxN,
        'mSF_median': np.median(mSF_list),
        'rtm_fraction': np.mean(rtm_list)}
    for i, val in enumerate(success):
        result[f'success{Nlist[i]}'] = val
    for i, val in enumerate(successR):
        result[f'successR{Nlist[i]}'] = val
    return result


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

def test3_N_for_efficacy_par(sensitivity,FAR,sz_min,sz_monthly_min,
                sz_longest_free,num_trials,num_patients,drug_effect,
                Nlist,use_baseline,correct_FAR):
    # note the eligibility rules come from Krauss et al 2022 cenobamate trial
    _require_baseline_eligibility(use_baseline)
    
    
    results = []
    results = Parallel(n_jobs=25)(delayed(process_efficacy)
                (sensitivity=sensitivity,FAR=FAR,num_patients=num_patients,correct_FAR=correct_FAR,
                use_baseline=use_baseline,sz_min=sz_min,sz_monthly_min=sz_monthly_min,
                sz_longest_free=sz_longest_free,drug_effect=drug_effect,Nlist=Nlist) 
                for num_trials in range(num_trials))
    
        # find the closest N value that produces success of 0.90
#    success = np.mean([[r[f'success{N}'] for N in Nlist] for r in results], axis=0)
    success = np.array([np.mean([r[f'success{N}'] for r in results]) for N in Nlist])
    successR = np.array([np.mean([r[f'successR{N}'] for r in results]) for N in Nlist])

    fraction_eligible = np.mean([r['fraction_eligible'] for r in results])

    target_success = 0.90
    closest_N_idx = np.argmin(np.abs(success - target_success))
    closest_N_for_90 = Nlist[closest_N_idx] if len(Nlist) > 0 else np.nan
    closest_N_idxR = np.argmin(np.abs(successR - target_success))
    closest_N_for_90R = Nlist[closest_N_idxR] if len(Nlist) > 0 else np.nan 

    return pd.DataFrame(results),closest_N_for_90,closest_N_for_90R,fraction_eligible


def test4_N_for_Type1error_par(sensitivity,FAR,sz_min,sz_monthly_min,
                sz_longest_free,num_trials,num_patients,
                Nlist,use_baseline,correct_FAR):
    _require_baseline_eligibility(use_baseline)
    
    # we assume no drug effect for type 1 error
    drug_effect = 0.0
    
    results = []
    results = Parallel(n_jobs=25)(delayed(process_efficacy)
                (sensitivity=sensitivity,FAR=FAR,num_patients=num_patients,correct_FAR=correct_FAR,
                use_baseline=use_baseline,sz_min=sz_min,sz_monthly_min=sz_monthly_min,
                sz_longest_free=sz_longest_free,drug_effect=drug_effect,Nlist=Nlist) 
                for num_trials in range(num_trials))
    
        # find the closest N value that produces success of 0.90
#    success = np.mean([[r[f'success{N}'] for N in Nlist] for r in results], axis=0)
    success = np.array([np.mean([r[f'success{N}'] for r in results]) for N in Nlist])
    successR = np.array([np.mean([r[f'successR{N}'] for r in results]) for N in Nlist])

    return pd.DataFrame(results),success,successR
