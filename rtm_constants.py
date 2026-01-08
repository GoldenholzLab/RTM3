# constants needed for RTM3 project

ELIGIBILITY_MONTHS = 2
DAYS_PER_MONTH = 30
TRIAL_BASELINE_MONTHS = 2
TRIAL_TEST_MONTHS = 3

EPSILON = 1e-10  # small constant to avoid divide by zero errors
COHORT_SIZE = 300*2# default cohort size

# for simple simulation
MONTHLY_STD = 2  # standard deviation for daily seizure counts in simplified model

# for plotting
HEATMAP_BINS = 30  # number of bins for heatmap plots

OUTPUT_FILENAME = 'output.csv'  # default output filename