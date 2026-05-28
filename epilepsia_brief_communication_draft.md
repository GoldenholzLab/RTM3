# Device sensitivity and false alarms can reshape regression to the mean in simulated epilepsy trials

**Article type:** Brief Communication  
**Running title:** Device errors and RTM  
**Authors:** Daniel M. Goldenholz<sup>1,2</sup>, Rohan Bhansali<sup>1,2</sup>, Brandon Westover<sup>2,3</sup> [affiliations to confirm]  
**Affiliations:**  
1. Department of Neurology, Beth Israel Deaconess Medical Center, Boston, Massachusetts, USA  
2. Harvard Medical School, Boston, Massachusetts, USA  
3. Department of Neurology, Massachusetts General Hospital, Boston, Massachusetts, USA  

**Corresponding author:** Daniel M. Goldenholz, Department of Neurology, Beth Israel Deaconess Medical Center, Boston, Massachusetts, USA. Address, telephone, fax, and email: [confirm before submission].  
**Key words:** seizure detection; clinical trials; regression to the mean; placebo response; seizure diary; false alarm rate  
**Manuscript counts:** Summary 180 words; main text 1,299 words; references 14; figures 2; tables 0; text pages [confirm after Word formatting].  

## Summary

Automated seizure detection devices are increasingly plausible tools for epilepsy trials, but no device is perfect. We used CHOCOLATES, a realistic seizure diary simulator, to examine how device sensitivity and false alarm rate (FAR) affect regression to the mean (RTM) and placebo median percentage change (MPC) in a simulated randomized trial design. For each device condition, 100,000 potential participants were generated; eligibility was assessed only during a 2-month baseline, followed by a 3-month test period. With FAR fixed at 0, reducing sensitivity from 100% to 10% increased the fraction of eligible participants exhibiting RTM from 38.2% to 64.8% and increased placebo MPC from 14.7% to 48.1%. With sensitivity fixed at 100% and expected FAR subtracted, increasing FAR from 0 to 1 alarm/day increased RTM from 38.2% to 53.2% and placebo MPC from 14.7% to 31.3%. Imperfect seizure detection can therefore change the apparent placebo response expected from trial eligibility criteria alone. Device operating characteristics should be incorporated into epilepsy trial simulations before objective seizure counts are used as trial endpoints.

## Short summary for table of contents

In simulated epilepsy trials, imperfect seizure detection altered regression to the mean and placebo median percentage change. Trial planning should model detector sensitivity and false alarm rate before device-derived seizure counts are used as endpoints.

## Introduction

Automated seizure detection is moving from an aspirational technology toward practical clinical and research use. Current devices are evaluated primarily by sensitivity, false alarm rate (FAR), detection latency, and usability, and recent guideline and meta-analytic work emphasizes that performance varies across seizure types and recording settings.<sup>1,2</sup> These metrics are usually interpreted at the level of patient safety, caregiver notification, or diary accuracy. Their consequences for randomized controlled trial (RCT) design are less developed.

Epilepsy trials are unusually vulnerable to apparent placebo response because seizure frequency is noisy, eligibility often selects patients during temporarily high seizure periods, and outcomes are defined as change from baseline.<sup>3-7</sup> Prior work has shown that regression to the mean (RTM) can be manipulated by eligibility criteria and can materially affect placebo response and trial power.<sup>3</sup> A separate line of work suggests that natural seizure-frequency variability can reproduce many features of placebo response even without a psychobiological placebo effect.<sup>5-7</sup>

Seizure detection devices change the measured seizure time series before trial eligibility and outcomes are calculated. Lower sensitivity can miss true seizures; FAR can add false events. Either process could change who qualifies for a trial and the apparent change from baseline. We asked how sensitivity and FAR alter RTM and placebo MPC in a simulated epilepsy RCT, using only baseline-period eligibility, a design choice aligned with current trial practice and the present analysis scope.

## Methods

We performed a simulation study using CHOCOLATES, an open-source seizure diary simulator that recapitulates empirical seizure diary properties, including heterogeneous long-term seizure rates, the relationship between mean seizure frequency and variability, seizure clustering, cycles, and interseizure constraints.<sup>8</sup> Analyses were performed in the RTM3 repository using `rtm_tools.py`, `plotting_tools.py`, and `for_rtm3_paper.ipynb`.

Each simulated participant had a 2-month baseline and 3-month test period. There was no active treatment effect. Eligibility was tested only during baseline. The baseline criteria were modeled after focal seizure trial eligibility rules used in cenobamate studies: mean baseline frequency at least 4 seizures/month, each baseline month at least 3 seizures, and no seizure-free interval longer than 25 days.<sup>9,10</sup>

Device sensitivity was implemented probabilistically at the daily seizure-count level: each true seizure had probability equal to the specified sensitivity of being detected. FAR was implemented as Poisson-distributed false alarms added to each day. For the primary displayed FAR analysis, expected false alarms were subtracted from monthly counts and from daily counts used for the seizure-free interval rule. This represents an idealized setting in which device FAR is known and is explicitly accounted for when deriving trial seizure counts.

Two device-performance sweeps were conducted. First, FAR was fixed at 0 alarms/day while sensitivity varied from 10% to 100% in 10% increments. Second, sensitivity was fixed at 100% while FAR varied from 0 to 1 alarm/day in 0.1 alarm/day increments. For each condition, 100,000 potential participants were generated. Among eligible participants, RTM was defined as a baseline seizure frequency above the participant's effective long-term seizure frequency and closer proximity of the test-period frequency to that long-term frequency. Placebo response was quantified as MPC, defined as 100 x (1 - test-period seizure frequency / baseline seizure frequency). Confidence intervals in Figure 1 are nonparametric 95% confidence intervals for the estimated median MPC, calculated from binomial order statistics around the sample median. They are not 95% intervals of individual patient percentage changes. No hypothesis tests were performed, because the analysis was a descriptive simulation sweep rather than a comparison of randomized groups. No missing-data procedure was required because all simulated diaries were complete.

This was a simulation study using synthetic participants only. No human participants, identifiable private information, or animals were used.

## Results

When FAR was 0, lower sensitivity increased both RTM and placebo MPC (Figure 1A,C). At 100% sensitivity, 33,435 of 100,000 potential participants were eligible, 38.2% of eligible participants met the RTM definition, and median placebo MPC was 14.7% (95% CI 13.9%-15.2%). At 50% sensitivity, 15,724 participants were eligible, RTM was 46.4%, and median placebo MPC was 25.9% (95% CI 25.0%-26.7%). At 10% sensitivity, only 1,660 participants were eligible, but RTM increased to 64.8%, and median placebo MPC increased to 48.1% (95% CI 45.5%-50.0%).

When sensitivity was 100%, increasing FAR also increased RTM and placebo MPC despite correction for expected false alarms (Figure 1B,D). At FAR 0.1 alarm/day, 35,377 participants were eligible, RTM was 40.7%, and median placebo MPC was 17.9% (95% CI 17.1%-18.5%). At FAR 0.5 alarms/day, 35,522 participants were eligible, RTM was 48.7%, and median placebo MPC was 26.8% (95% CI 26.2%-27.5%). At FAR 1 alarm/day, 35,500 participants were eligible, RTM was 53.2%, and median placebo MPC was 31.3% (95% CI 30.3%-33.3%).

Thus, both missed seizures and false alarms changed the placebo response expected from natural variability and baseline selection alone. The direction was consistent in the primary displayed analysis: lower sensitivity and higher FAR increased the RTM fraction and increased placebo MPC.

## Discussion

This simulation suggests that device performance characteristics should be treated as trial-design parameters, not merely as validation metrics. In a conventional baseline-versus-treatment design, the measured baseline seizure rate determines eligibility and also serves as the denominator for percent-change outcomes. Imperfect detection changes both. Lower sensitivity enriched the eligible cohort for unusually high measured baseline periods among people whose true long-term seizure frequency was lower, producing larger subsequent apparent improvement. FAR had a similar effect even when expected false alarms were subtracted, because false alarms could still alter daily and monthly eligibility behavior and trial selection.

These findings extend prior RTM work in epilepsy RCTs.<sup>3</sup> The earlier analysis showed that eligibility criteria can increase or decrease RTM and placebo response. The present analysis adds that the measurement system itself can alter the same pathway. If future trials use device-derived seizure counts, then two trials with identical biological treatment effects and identical eligibility rules could have different apparent placebo MPC solely because the devices, algorithms, or postprocessing rules differ.

The practical implication is not that device-derived endpoints should be avoided. Objective detection could reduce under-reporting and may improve statistical power when more complete seizure counts are available.<sup>11</sup> Rather, device operating characteristics should be incorporated into trial simulations before design finalization. Sponsors and investigators should prespecify how sensitivity, FAR, and FAR correction will enter eligibility, baseline counts, outcome counts, and missing-data rules. Device validation should also report performance in the seizure types and clinical settings that match the planned trial population.<sup>1,2</sup>

Several limitations matter. First, CHOCOLATES is a simulator, not a biological model of every epilepsy syndrome. Its value is that it reproduces multiple statistical properties of seizure diaries and trial-like placebo responses, but it remains an abstraction.<sup>8</sup> Second, sensitivity and FAR were modeled simply and uniformly across participants and time. Real detectors may have patient-specific performance, seizure-type-specific sensitivity, state dependence, clustering of false alarms, and periods of device nonwear or signal deficiency. Third, the FAR correction was idealized. In practice, false alarms may not be known at the individual-event level, and correction rules could affect eligibility and outcomes differently. Fourth, the confidence intervals in Figure 1 quantify uncertainty in the estimated median MPC from simulation, not the dispersion of patient-level responses.

Seizure detection devices may make epilepsy trials more objective. They will also make trial statistics depend on algorithm performance. RTM and placebo response should therefore be evaluated jointly with sensitivity and FAR whenever device-derived seizure counts are used for epilepsy RCT eligibility or endpoints.

## Acknowledgments and funding

Daniel M. Goldenholz (DMG) was funded by the National Institutes of Health (NIH) K23NS124656, 1R21NS142800 and the American Board of Psychiatry and Neurology (ABPN).

## Disclosure of conflicts of interest

DMG has been provided speaker fees from Harvard Medical School, the American Academy of Neurology, the American Epilepsy Society, the American Clinical Neurophysiology Society, the National Neurotrauma Society, AI in Epilepsy and Neurology, Florida Epilepsy Alliance, and the University of Texas at Austin. He also previously has been a paid consultant for Neuro Event Labs, IDR, LivaNova, Health Advances, Duke University, Bloom Insights, and Wiley. He has received grants from NIH, ABPN, Beth Israel Deaconess Medical Center, and the Lions Club. Rohan Bhansali and Brandon Westover have no conflicts of interest to disclose. [Author disclosures should be confirmed before submission.]

## Ethical publication statement

We confirm that we have read the Journal's position on issues involved in ethical publication and affirm that this report is consistent with those guidelines.

## Data and code availability

The simulations were performed using code in the RTM3 repository. The exact plotted results are archived in `rtm_test123_mpc_ci_results.csv`; the primary notebook is `for_rtm3_paper.ipynb`; the figure-generation code is in `plotting_tools.py`; and the simulation functions are in `rtm_tools.py`. CHOCOLATES is available from the GoldenholzLab GitHub repository as described previously.<sup>8</sup>

## Author contributions

Daniel M. Goldenholz conceived the analysis, wrote and revised analysis code, interpreted results, and drafted the manuscript. Rohan Bhansali contributed to interpretation and manuscript revision. Brandon Westover contributed to interpretation, simulation framing, and manuscript revision. [Confirm and revise author contribution language before submission.]

## References

1. Beniczky S, Wiebe S, Jeppesen J, et al. Automated seizure detection using wearable devices: A clinical practice guideline of the International League Against Epilepsy and the International Federation of Clinical Neurophysiology. Clin Neurophysiol 2021;132:1173-1184.
2. Naganur V, Sivathamboo S, Chen Z, et al. Automated seizure detection with noninvasive wearable devices: A systematic review and meta-analysis. Epilepsia 2022;63:1930-1941.
3. Goldenholz DM, Goldenholz EB, Kaptchuk TJ. Quantifying and controlling the impact of regression to the mean on randomized controlled trials in epilepsy. Epilepsia 2023;64:2635-2643.
4. Goldenholz DM, Goldenholz SR. Response to placebo in clinical epilepsy trials: Old ideas and new insights. Epilepsy Res 2016;122:15-25.
5. Goldenholz DM, Moss R, Scott J, et al. Confusing placebo effect with natural history in epilepsy: A big data approach. Ann Neurol 2015;78:329-336.
6. Goldenholz DM, Strashny A, Cook M, et al. A multi-dataset time-reversal approach to clinical trial placebo response and the relationship to natural variability in epilepsy. Seizure 2017;53:31-36.
7. Romero J, Larimer P, Chang B, et al. Natural variability in seizure frequency: Implications for trials and placebo. Epilepsy Res 2020;162:106306.
8. Goldenholz DM, Westover MB. Flexible realistic simulation of seizure occurrence recapitulating statistical properties of seizure diaries. Epilepsia 2023;64:396-405.
9. Chung SS, French JA, Kowalski J, et al. Randomized phase 2 study of adjunctive cenobamate in patients with uncontrolled focal seizures. Neurology 2020;94:e2311-e2322.
10. Krauss GL, Klein P, Brandt C, et al. Safety and efficacy of adjunctive cenobamate (YKP3089) in patients with uncontrolled focal seizures: A multicentre, double-blind, randomised, placebo-controlled, dose-response trial. Lancet Neurol 2020;19:38-48.
11. Goldenholz DM, Tharayil JJ, Kuzniecky R, et al. Simulating clinical trials with and without intracranial EEG data. Epilepsia Open 2017;2:156-161.
12. Rheims S, Perucca E, Cucherat M, et al. Factors determining response to antiepileptic drugs in randomized controlled trials. A systematic review and meta-analysis. Epilepsia 2011;52:219-233.
13. Zaccara G, Giovannelli F, Schmidt D. Placebo and nocebo responses in drug trials of epilepsy. Epilepsy Behav 2015;43:128-134.
14. McDonald CJ, Mazzuca SA, McCabe GP Jr. How much of the placebo "effect" is really statistical regression? Stat Med 1983;2:417-427.

## Figure legends

**Figure 1. Device sensitivity, false alarm rate, regression to the mean, and placebo median percentage change.** File: `Goldenholz-fig1.tif` (source PNG: `sensitivity_and_FAR_vs_RTM_mpc_ci.png`). Simulations used 100,000 potential participants per condition and baseline-only eligibility. Panels A and C vary sensitivity with FAR fixed at 0 alarms/day. Panels B and D vary FAR with sensitivity fixed at 100%, with expected false alarms subtracted. RTM, regression to the mean; FAR, false alarm rate; MPC, median percentage change. Error bars show 95% confidence intervals for the estimated median MPC.

**Figure 2. Conceptual mechanism linking detector error to apparent placebo response.** File: `Goldenholz-fig2.tif` (source PNG: `rtm_far_sensitivity_explainer.png`). Lower sensitivity and higher FAR change the measured seizure-count series used for trial eligibility and baseline comparison. Because eligibility selects higher measured baseline seizure counts, detector error can enrich for participants whose later measured counts regress toward their long-term seizure frequency, increasing apparent placebo response.
