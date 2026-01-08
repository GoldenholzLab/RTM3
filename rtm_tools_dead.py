


# write a generator for simplified patients
def simplified_patient_generator(mSF=None, scale=CONST.MONTHLY_STD):
    if mSF is None:
        mSF = choose_seizure_frequency()

    # normally distrubuted daily seizure counts
    diaryPRE = np.random.normal(loc=mSF, scale=scale, size=1)
    diaryBASELINE = np.random.normal(loc=mSF, scale=scale, size=1)
    diaryTEST = np.random.normal(loc=mSF, scale=scale, size=1)
    
    # make sure these numbers are non-negative
    diaryPRE[diaryPRE < 0] = 0
    diaryBASELINE[diaryBASELINE < 0] = CONST.EPSILON
    diaryTEST[diaryTEST < 0] = 0

    return diaryPRE, diaryBASELINE, diaryTEST, mSF

# generate a cohort of simplified patients
def generate_simplified_cohort(num_patients=CONST.COHORT_SIZE, STD=CONST.MONTHLY_STD, eligibility_min=0):
    cohort_data = []
    per_arm = num_patients // 2

    for _ in range(per_arm):
        diaryPRE, diaryBASELINE, diaryTEST, mSF = simplified_patient_generator(mSF=None, scale=STD)
        # note I don't treat this arm
        cohort_data.append({
            'diaryPRE': diaryPRE[0],
            'diaryBASELINE': diaryBASELINE[0],
            'diaryTEST': diaryTEST[0],
            'mSF': mSF[0],
            'PC': 100*(1 - diaryTEST[0] / diaryBASELINE[0])
        })
    
    cohort_data2 = keepEligible(cohort_data, eligibility_min, use_baseline=True)

    return pd.DataFrame(cohort_data2)

def keepEligible(cohort_data, eligibility_min, use_baseline=True):
    eligible_data = []
    for patient in cohort_data:
        if use_baseline:
            freq = patient['diaryBASELINE']
        else:
            freq = patient['diaryPRE']
        if freq >= eligibility_min:
            eligible_data.append(patient)
    return eligible_data


def generated_one_complex_patient(use_baseline=True, sensitivity=1.0, FAR=0.0):
    if use_baseline:
        preM = CONST.ELIGIBILITY_MONTHS - CONST.TRIAL_BASELINE_MONTHS
    else:
        preM = CONST.ELIGIBILITY_MONTHS


    daysNeeded = (preM + CONST.TRIAL_BASELINE_MONTHS + CONST.TRIAL_TEST_MONTHS) * CONST.DAYS_PER_MONTH
    mSF = choose_seizure_frequency()
    dailyDiary= simulator_base(sampRATE=1, number_of_days=daysNeeded, defaultSeizureFreq=mSF)
    monthlyDiary = downsample(x=dailyDiary,byHowmuch=CONST.DAYS_PER_MONTH)

    # address sensitivity and FAR
    monthlyDiary = (sensitivity * monthlyDiary).astype(int)
    monthlyDiary = monthlyDiary + np.random.poisson(FAR*CONST.DAYS_PER_MONTH, size=len(monthlyDiary))
    # apply correction factor
    monthlyDiary = monthlyDiary - int(FAR*CONST.DAYS_PER_MONTH)
    monthlyDiary[monthlyDiary < 0] = 0  # ensure no negative seizure


    diaryPre = np.mean(monthlyDiary[:preM])
    diaryBaseline = np.mean(monthlyDiary[preM:(preM + CONST.TRIAL_BASELINE_MONTHS)])
    if diaryBaseline < CONST.EPSILON:
        diaryBaseline = CONST.EPSILON  # avoid divide by zero
    diaryTest = np.mean(monthlyDiary[(preM + CONST.TRIAL_BASELINE_MONTHS):])    
    return diaryPre, diaryBaseline, diaryTest, mSF

# generate complex cohort
def generate_complex_cohort(num_patients=CONST.COHORT_SIZE,eligibility_min=0,use_baseline=True,
                            sensitivity=1.0,FAR=0.0):

    cohort_data = []
    per_arm = num_patients // 2
    for _ in range(per_arm):
        diaryPre, diaryBASELINE, diaryTEST, mSF = generated_one_complex_patient(use_baseline=use_baseline, sensitivity=sensitivity, FAR=FAR)
        cohort_data.append({
                    'diaryPRE': diaryPre,
                    'diaryBASELINE': diaryBASELINE,
                    'diaryTEST': diaryTEST,
                    'mSF': mSF[0],
                    'PC': 100*(1 - diaryTEST / diaryBASELINE)
                })
    cohort_data2 = keepEligible(cohort_data, eligibility_min, use_baseline=use_baseline)
    return pd.DataFrame(cohort_data2)