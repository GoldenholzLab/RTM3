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

# modify the test period seizure rate based on treatment effect
def treat_one_patient(diaryTEST, treatment_effect):
    # inputs
    #   diaryTEST - numpy array of monthly seizure counts
    #   treatment_effect - reduction in monthly seizure frequency due to treatment
    # outputs
    #   treated_diaryTEST - numpy array of monthly seizure counts after treatment

    treated_diaryTEST = diaryTEST * (1 - treatment_effect)
    treated_diaryTEST[treated_diaryTEST < 0] = 0  # ensure no negative seizure counts
    return treated_diaryTEST

# wrapper for CHCOCOLATES
def simulate_one(daysNeeded)
    mSF = get_mSF(requested_msf=-1)
    dailyDiary= simulator_base(sampRATE=1, number_of_days=daysNeeded, defaultSeizureFreq=mSF)
    return dailyDiary, mSF

def apply_sens_and_far_to_dairy(dailyDiary, sensitivity, FAR):
    # inputs
    #   dailyDiary - numpy array of daily seizure counts
    #   sensitivity - sensitivity of the detection algorithm
    #   FAR - false alarm rate (false positives per day)
    # outputs
    #   detectedDiary - numpy array of daily detected seizure counts
        # apply sensitivity probabilistically - each seizure has probability 'sensitivity' of being detected
    if sensitivity < 1.0:
        dailyDiary = np.array([np.random.binomial(n=int(count), p=sensitivity) if count > 0 else 0 for count in dailyDiary])
    # apply FAR probabilistically    
    if FAR > 0.0:
        dailyDiary += np.random.poisson(FAR,size=len(dailyDiary))

    return dailyDiary

def nake_a_quick_test_of_apply_sens_and_far_to_dairy():
    # quick test of the apply_sens_and_far_to_dairy function
    dailyDiary = np.ones(1000)*10
    correct_FAR = True
    for sensitivity in np.arange(0.0,1.1,0.1):
        for FAR in [0.0, 0.5, 1.0]:
            detectedDiary = apply_sens_and_far_to_dairy(dailyDiary, sensitivity, FAR, correct_FAR)
            est_FAR = np.mean(detectedDiary) - np.mean(dailyDiary)

    sensitivity = 0.5
    FAR = 1.0
    detectedDiary = apply_sens_and_far_to_dairy(dailyDiary, sensitivity, FAR)
    print("Original diary: ", dailyDiary)
    print("Detected diary: ", detectedDiary)