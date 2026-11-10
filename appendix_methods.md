# RTM3 scientific appendix

This appendix reports every simulated condition, selection diagnostics, correction comparisons, and the fixed-cohort sensitivity analysis. These additional analyses were performed during manuscript revision, not prespecified. All counts are simulated.

RTM denotes regression to the mean. False-alarm rate (FAR) is the expected number of false alarms per day.

The analyses ask how detector performance and the enrollment rule change apparent improvement when there is no treatment effect. Observed RTM is the proportion meeting a strict trajectory criterion. Median percentage change (MPC) measures the median percentage reduction from baseline. These are different outcomes and can move in different directions.

**The fixed perfect-detector cohort is hypothetical and completely unrealistic as an enrollment strategy for a trial whose eligibility is determined by the imperfect detector being evaluated.**

Some retained participants would have failed that detector's entry criteria. The comparison tests how the result depends on changing enrollment. This qualification concerns this particular enrollment strategy; it does not imply that every fixed-cohort trial design is unrealistic or that the diagnostic comparison is mathematically invalid.

## Shared participant characteristics and reproducibility

Each of 100,000 potential participants has a 60-day baseline, a 90-day test period, and no treatment effect. A separate 1,080-day realization estimates effective long-term frequency. A simulated month is 30 days throughout. The two realizations share the sampled monthly rate, cluster status and settings, and cycle frequencies and amplitudes. Their random count realizations differ. The long-term diary is not a continuation of the trial diary, and shared cycle parameters do not imply a shared daily intensity trajectory.

At one sample per day, the unchanged CHOCOLATES source cannot represent 12-hour or 24-hour within-day cycles. Its configured weekly component has zero amplitude. Longer cycles, when drawn, remain possible. The daily study therefore does not reproduce the simulator's full within-day rhythmic behavior. Empirical diaries used to develop CHOCOLATES can contain reporting error; treating generated events as latent biological truth is a modeling assumption.

The base seed is 20260813, the nested detection seed is 20260814, and the nested false-alarm seed is 20260815. Lower-sensitivity detections are contained in higher-sensitivity realizations. Higher false-alarm conditions add independent Poisson increments. Every eligibility/correction comparison uses the same latent participants and the same detector realization at that condition. Intervals quantify Monte Carlo sampling uncertainty conditional on this model, not uncertainty about its assumptions or real devices.

## Observations and monthly correction

For participant $i$ and day $d$, let $S_{id}$ be latent seizures, $s$ sensitivity, and $\lambda$ the false-alarm rate in alarms/day:

$$D_{id}\mid S_{id}\sim\operatorname{Binomial}(S_{id},s),\quad F_{id}\sim\operatorname{Poisson}(\lambda),\quad Y_{id}=D_{id}+F_{id}.$$

False alarms are independent of latent counts. Let $Y_{im}$ be the observed total in a 30-day month. Expected-alarm correction gives

$$C_{im}=\max(0,Y_{im}-30\lambda).$$

Here $D_{id}$ counts detected latent seizures, $F_{id}$ counts false alarms, and $Y_{id}$ counts all detector events. The corrected monthly total $C_{im}$ subtracts the expected number of alarms from observed monthly total $Y_{im}$ and floors negative values at zero.

Before flooring, subtracting the known expectation leaves mean-zero false-alarm error in the unselected population. It does not identify individual alarms. Flooring and selection can change the residual error among selected participants. This correction targets sensitivity-scaled latent counts, but it is not an unbiased estimator after flooring.

## Eligibility and additive thresholds

The main rule uses only the two baseline months:

$$E_i=\mathbf1\{C_{i1}\ge3,\ C_{i2}\ge3,\ (C_{i1}+C_{i2})/2\ge4\}.$$

The additive alternative selects on raw observed counts:

$$E_i^{\mathrm{add}}=\mathbf1\{Y_{i1}\ge3+30\lambda,\ Y_{i2}\ge3+30\lambda,\ (Y_{i1}+Y_{i2})/2\ge4+30\lambda\}.$$

The indicator $\mathbf1\{\cdot\}$ equals one when every condition inside the braces holds and zero otherwise. Thus $E_i=1$ means that participant $i$ is eligible; $E_i^{\mathrm{add}}$ is the same decision expressed on the raw-count scale. All baseline inequalities include equality.

The equivalent two-month total threshold is $8+60\lambda$, retaining the separate monthly minimum. Because each eligible corrected month is at least three, neither is on the zero floor. Subtracting $30\lambda$ from each raw month therefore yields exactly the primary inequalities. Thus $E_i^{\mathrm{add}}=E_i$ for every participant. This is an algebraic identity, not only a simulation result. The notebook verifies mask and outcome equality at all 11 false-alarm rates. Both analyses retain corrected outcome counts. Raising raw thresholds changes how the rule is expressed but does not improve precision or remove random false-alarm error.

## Effective long-term reference and strict endpoint

Let $L_i$ be the integer seizure total in the separate 36-month realization:

$$f_i^*=L_i/36,\qquad r_i=s f_i^*.$$

The effective long-term frequency differs from the sampled rate parameter because the realization includes count variability, clustering, and configured cycles. For mean monthly measured baseline $B_i$ and test frequency $T_i$, observed RTM is

$$R_i=\mathbf1\{B_i>r_i,\quad B_i-r_i>|T_i-r_i|\}.$$

Here $f_i^*$ is effective latent seizures per month and $r_i$ is its sensitivity-scaled reference. The measured frequencies $B_i$ and $T_i$ are the averages over two baseline months and three test months, respectively. The indicator $R_i$ requires baseline to exceed the reference and test to be strictly closer to it. Its prevalence is the number with $R_i=1$ divided by all selected participants.

Both comparisons are strict. Equal distances do not count. The canonical grid uses an exact common integer scale for the two-month baseline, three-month test, rational sensitivity, and 36-month reference. Retaining $L_i$ prevents floating-point errors at true ties. Genuinely noninteger inputs outside that grid use a documented relative equality tolerance of $10^{-10}$.

The separate uncorrected-outcome comparison uses its observation-scale reference:

$$r_i^{\mathrm{raw}}=s f_i^*+30\lambda.$$

That comparison uses raw counts for eligibility and outcomes with the original fixed count thresholds. It differs from the additive-threshold analysis, which retains corrected outcomes. Comparing raw false-alarm-contaminated outcomes with $s f_i^*$ would mix count scales and is not done here. Exact oracle removal subtracts realized labeled false alarms and uses the latent-count reference. It is invariant across false-alarm conditions by construction.

The primary aggregate observed-RTM analysis describes trajectories. It does not identify an isolated participant-level RTM Type 3 component or divide overall RTM into additive causal parts. The additional paired analysis below compares measured and latent classifications within each detector-selected cohort and reports an operational measurement-only category, with separate reverse discordance.

## Percentage change and denominators

For positive measured baseline:

$$PC_i=100(1-T_i/B_i),\qquad\mathrm{MPC}=\operatorname{median}_{i:E_i=1,\ B_i>0}(PC_i).$$

Here $PC_i$ is a participant's percentage reduction from measured baseline, so positive values mean fewer measured seizures during test. MPC is the median of those defined reductions; it is not the percentage of participants showing observed RTM. For fixed membership, the selection indicator in this median refers to perfect-detector enrollment, while $B_i$ and $T_i$ still come from the detector condition being evaluated.

All main-analysis participants have positive baselines. Fixed-cohort participants can have zero measured baseline after detector performance changes. Their percentage change is undefined, regardless of test count. Tables report the full selected count, zero-baseline count, and positive-baseline MPC denominator. No pseudocount is added. Observed RTM always uses the full selected denominator, including zero-baseline participants, who do not meet its criterion.

MPC intervals are nonparametric binomial-order-statistic 95% intervals. CSVs also contain Wilson intervals for selection and observed RTM. No hypothesis tests or between-condition confidence intervals are reported. Conditions share participants, so pointwise intervals are not independent-comparison intervals.

## Fixed versus reselected cohorts

Reselection applies eligibility to the condition-specific detector baseline counts. It models that detector determining enrollment. The fixed cohort consists of the 35,393 patients who met eligibility under perfect detection, $s=1$ and $\lambda=0$. Their identities are retained regardless of whether they would qualify under an imperfect detector. Both strategies share the latent seizure histories and the exact same detector realizations at each condition. Membership is fixed, not successful eligibility under each detector.

This is a hypothetical diagnostic comparison of dependence on changing enrollment. It is not a realistic prediction of enrollment in a detector-selected trial, an additive causal decomposition, or an isolated RTM Type 3 estimate. At 10% sensitivity, 33,213 retained participants fail the current detector's criteria; at one false alarm/day, 11,351 fail them. These failures are distinct from the 6,182 and 538 zero measured baselines excluded from MPC, respectively. Observed RTM retains all 35,393 patients in its denominator.

Fixed-cohort outcomes describe changes in measured trajectories without changing initial enrollment. However, MPC remains conditional on positive measured baseline, and that subset can vary. The selected count is constant, but its MPC denominator need not be. This limits interpretation of low-sensitivity fixed-cohort MPC.

## The 25-day rule

Let $L_i(Z)$ be the longest consecutive zero-day run in the 60-day baseline. The interval sensitivity adds $L_i(Z)\le25$ to monthly criteria. The original daily correction is

$$Z_{id}=\max(0,Y_{id}-\operatorname{round}(\lambda)).$$

Here $Z_{id}$ is the daily count used only to determine zero-day runs; $L_i(Z)$ is the longest such run. The maximum run is inclusive at 25 days. Monthly enrollment criteria and outcome counts continue to use $C_{im}$.

NumPy nearest-even rounding subtracts zero from zero through 0.5/day and one from 0.6 through 1.0/day. This discontinuity changes zero-day classification without identifying false alarms. Alternatives use observed days, $Z_{id}=Y_{id}$, or oracle detected-seizure days, $Z_{id}=D_{id}$. Only the interval component changes; monthly outcomes retain expected-alarm correction. Oracle intervals and exact oracle alarm removal are separate diagnostics.

The 30-day months adapt clinical criteria rather than reproducing four-week trial windows exactly. Dropping the interval rule removes the need to assign corrected events to days, but it does not eliminate monthly noise or flooring.
