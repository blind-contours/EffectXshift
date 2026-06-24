
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`EffectXshift`

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/EffectXshift/actions/workflows/r.yml/badge.svg?branch=main)](https://github.com/blind-contours/EffectXshift/actions/workflows/r.yml)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/EffectXshift/master.svg)](https://codecov.io/github/blind-contours/EffectXshift?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447) -->
<!-- badges: end -->

> Effect Modification Identification and Estimation using Data-Adaptive
> Stochastic Interventions **Authors:** [David
> McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What’s `EffectXshift`?

The `EffectXshift` R package offers an approach which identifies and
estimates the differential impact of intervention on a mixed exposure on
an outcome. We define effect modification as a contrast of
stochastic-shift effects: the expected change in the outcome under an
additive shift of one exposure in a subregion of the covariate space
compared to the complementary region of that space. Stochastic
interventions or exposure changes depend on naturally observed values,
as described in past literature (Dı́az and van der Laan 2012; Haneuse
and Rotnitzky 2013).

Our target parameter is:

$$
\begin{align*}
(\hat A_i, \hat V)
&= \underset{A_i \in \boldsymbol{A}, \, V \subseteq W}{\arg\max}
\left\{
\Psi_i(V) - \Psi_i(V^c)
\right\}, \\
\Psi_i(V)
&= E\left[
E\{Y \mid A_i + \delta_i, A_{\setminus i}, W\}
- E\{Y \mid A_i, A_{\setminus i}, W\}
\mid W \in V
\right].
\end{align*}
$$

This finds the exposure-covariate region combination where the effect of
an intervention is maximally different.

## Under the Hood

`EffectXshift` first identifies effect modification using g-computation
through a Super Learner (ensemble of machine learning algorithms). We
get the exposure shift for each individual for each exposure. We then
regress these effects onto the covariate space using a custom decision
tree algorithm that aims to find the region in the covariate space where
the effect of shifting is maximally different. This algorithm therefore
selects the exposure from a mixture paired with a region in the
covariate space that maximizes the subpopulation intervention effect, or
rather, finds the exposure-covariate region combination where the effect
of intervention is maximally different.

## K-Fold Cross-Validation

The package uses a k-fold cross-validation framework for the
data-adaptive parameter. The process begins by partitioning the data into
parameter-generating and estimation samples. In the parameter-generating
sample, we identify our exposure-covariate region with maximum effect
modification using a machine learning g-computation framework. In the
estimation sample we then estimate the stochastic-shift effect in the
discovered region and its complement using the doubly robust estimator
TMLE, which allows us to construct confidence intervals for the held-out
estimates.

By using EffectXshift, users get access to a tool that offers both
k-fold specific and aggregated results for the maximal effect
modification, ensuring that researchers can glean the most information
from their data. For a more in-depth exploration, there’s an
accompanying vignette.

## Inputs

To utilize the package, users need to provide vectors for exposures,
covariates, and outcomes. They also specify the respective $\delta$ for
each exposure (indicating the degree of shift); unnamed delta vectors are
matched to exposure columns in order. The `top_n` parameter defines the
top number of effect modification exposure-covariate pairs. A detailed
guide is provided in the vignette. With these inputs, `EffectXshift`
processes the data and delivers tables showcasing fold-specific results
and aggregated outcomes, allowing users to glean insights effectively.

`EffectXshift` also incorporates features from the `sl3` package (Coyle,
Hejazi, Malenica, et al. 2022), facilitating ensemble machine learning
in the estimation process. If the user does not specify any stack
parameters, `EffectXshift` will automatically create an ensemble of
machine learning algorithms that strike a balance between flexibility
and computational efficiency.

## Where This Fits

`EffectXshift` is not intended to replace all heterogeneous-effect or
mixture methods. Its target is more specific: interpretable discovery and
held-out estimation of covariate regions where stochastic-shift effects
differ for one component of a continuous exposure mixture.

| Literature | Typical target | How `EffectXshift` differs |
|---|---|---|
| Stochastic interventions / modified treatment policies (Díaz and van der Laan 2012; Haneuse and Rotnitzky 2013) | Population mean under an exposure shift | Adds data-adaptive exposure and covariate-region discovery for effect modification |
| Honest causal trees and forests (Athey and Imbens 2016; Wager and Athey 2018) | Heterogeneous effects, often for binary or simple treatments | Targets additive stochastic shifts for continuous mixture components and reports interpretable selected regions |
| Meta-learners for CATE (Künzel et al. 2019) | Individualized conditional treatment effects | Uses fold-split discovery and TMLE/onestep estimation for a pre-specified stochastic-shift estimand |
| Exposure-mixture methods such as quantile g-computation (Keil et al. 2020) | Overall mixture effects or component weights | Searches for the exposure-region pair with maximal shift-effect contrast |

Because the exposure/region is selected from the data, users should
review the new selection and positivity diagnostics before interpreting
pooled estimates. Unstable fold-level selections suggest an exploratory
finding; large clever covariates suggest weak practical support for the
requested exposure shift.

------------------------------------------------------------------------

## Installation

*Note:* Because the `EffectXshift` package (currently) depends on `sl3`
that allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `EffectXshift` package is not
available on CRAN and must be downloaded here.

``` r
remotes::install_github("tlverse/sl3@master")
```

Make sure `sl3` installs correctly then install `EffectXshift`

``` r
remotes::install_github("blind-contours/EffectXshift@main")
```

## Trialist Guides

The pkgdown site includes trial-focused pages that are more practical
than the full methods vignette:

- [Trialist Quick
  Start](https://blind-contours.github.io/EffectXshift/articles/trialist-quick-start.html):
  estimand checklist, required inputs, output map, and reporting
  template.
- [Simulated Randomized Trial
  Walkthrough](https://blind-contours.github.io/EffectXshift/articles/simulated-rct-trialist.html):
  simulated trial data, `rct = TRUE` analysis, result tables,
  diagnostics, and interpretation.
- [Full Methods
  Vignette](https://blind-contours.github.io/EffectXshift/articles/EffectXshift-vignette.html):
  broader stochastic-shift and mixed-exposure workflow.

For a randomized trial, start with the quick-start page, then run
through the simulated RCT article before applying the method to a real
SAP or exploratory subgroup analysis.

## Example

To illustrate how `EffectXshift` may be used to ascertain the effect of
a mixed exposure, consider the following example:

``` r
library(EffectXshift)
library(devtools)
#> Loading required package: usethis
library(kableExtra)
library(sl3)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:kableExtra':
#> 
#>     group_rows
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union


seed <- 429153
set.seed(seed)
```

We will simulate some data with effect modification to see if we both
find the correct exposure-covariate combination that has the maximum
effect modification and compare our estimates in each level of the
modifier for an intervention on the exposure.

``` r
n <- 2000 # Number of observations

# Generate binary covariates
W3 <- rbinom(n, 1, 0.5) # 0 for one sex, 1 for another
W2 <- rbinom(n, 1, 0.5) # Additional binary covariate

# Generate confounders
W1 <- rbinom(n, 1, 0.5)

# Generate exposures influenced by confounders
A1 <- rnorm(n, mean = 0.5 * W3 + 0.3 * W2 + 0.4 * W1) # Continuous exposure with significant interaction with Sex
A2 <- rnorm(n, mean = 0.3 * W2 + 0.3 * W1) # Continuous exposure without significant interaction
A3 <- rnorm(n, mean = 0.2 * W1) # Continuous exposure without significant interaction

# Define effect sizes
beta_A1 <- 1 # Base effect of A1 on Y
beta_A2 <- 0.5 # Effect of A2 on Y
beta_A3 <- 0.2 # Effect of A3 on Y
beta_W3 <- 0.5 # Effect of Sex on Y
beta_W2 <- -0.3 # Effect of W2 on Y
beta_W1 <- 0.4 # Effect of W1 on Y
interaction_A1_W3 <- 2 # Interaction effect between A1 and Sex, ensuring maximal impact difference

# Simulate outcome Y
Y <- 2 + beta_A1 * A1 + beta_A2 * A2 + beta_A3 * A3 + beta_W3 * W3 + beta_W2 * W2 + beta_W1 * W1 +
  interaction_A1_W3 * A1 * W3 + rnorm(n) # Include noise

# Assume an intervention that reduces A1 by a fixed amount for the shift analysis
shift_amount <- -0.5
A1_shifted <- A1 + shift_amount

# Recalculate Y assuming the shift in A1
Y_shifted <- 2 + beta_A1 * A1_shifted + beta_A2 * A2 + beta_A3 * A3 + beta_W3 * W3 + beta_W2 * W2 + beta_W1 * W1 + interaction_A1_W3 * A1_shifted * W3 + rnorm(n)

# Data frame with shifted outcomes
data_shifted <- data.frame(W3, W2, W1, A1, A1_shifted, A2, A3, Y, Y_shifted)

effect_summary <- data_shifted %>%
  mutate(Difference = Y_shifted - Y) %>%
  group_by(W3, W2, W1) %>%
  summarise(
    Avg_Difference = mean(Difference),
    .groups = "drop"
  )

print(effect_summary)
#> # A tibble: 8 × 4
#>      W3    W2    W1 Avg_Difference
#>   <int> <int> <int>          <dbl>
#> 1     0     0     0         -0.549
#> 2     0     0     1         -0.388
#> 3     0     1     0         -0.362
#> 4     0     1     1         -0.606
#> 5     1     0     0         -1.43 
#> 6     1     0     1         -1.60 
#> 7     1     1     0         -1.56 
#> 8     1     1     1         -1.57

data_shifted$Y_diff <- data_shifted$Y_shifted - data_shifted$Y

# Aggregate differences by levels of Sex
avg_diff_W3 <- aggregate(Y_diff ~ W3, data = data_shifted, mean)

# Aggregate differences by levels of W2
avg_diff_W2 <- aggregate(Y_diff ~ W2, data = data_shifted, mean)

# Aggregate differences by levels of W1
avg_diff_W1 <- aggregate(Y_diff ~ W1, data = data_shifted, mean)

print(avg_diff_W3)
#>   W3     Y_diff
#> 1  0 -0.4745651
#> 2  1 -1.5385141
print(avg_diff_W2)
#>   W2     Y_diff
#> 1  0 -0.9881843
#> 2  1 -1.0236184
print(avg_diff_W1)
#>   W1     Y_diff
#> 1  0 -0.9980014
#> 2  1 -1.0147684
```

``` r
w <- data_shifted[, c("W1", "W2", "W3")]
a <- data_shifted[, c("A1", "A2", "A3")]
y <- data_shifted$Y


deltas <- list(
  "A1" = -0.5, "A2" = -0.5, "A3" = -0.5
)
head(data_shifted) %>%
  kbl(caption = "Sim Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Sim Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
W3
</th>
<th style="text-align:right;">
W2
</th>
<th style="text-align:right;">
W1
</th>
<th style="text-align:right;">
A1
</th>
<th style="text-align:right;">
A1_shifted
</th>
<th style="text-align:right;">
A2
</th>
<th style="text-align:right;">
A3
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
Y_shifted
</th>
<th style="text-align:right;">
Y_diff
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.6561885
</td>
<td style="text-align:right;">
0.1561885
</td>
<td style="text-align:right;">
1.2513750
</td>
<td style="text-align:right;">
-0.0302229
</td>
<td style="text-align:right;">
5.7034805
</td>
<td style="text-align:right;">
4.5714986
</td>
<td style="text-align:right;">
-1.1319819
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
-0.0663319
</td>
<td style="text-align:right;">
-0.5663319
</td>
<td style="text-align:right;">
-0.7459642
</td>
<td style="text-align:right;">
-0.2913230
</td>
<td style="text-align:right;">
1.1732674
</td>
<td style="text-align:right;">
0.1622237
</td>
<td style="text-align:right;">
-1.0110437
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.7725338
</td>
<td style="text-align:right;">
-1.2725338
</td>
<td style="text-align:right;">
-1.7225836
</td>
<td style="text-align:right;">
-0.4569841
</td>
<td style="text-align:right;">
-0.4378984
</td>
<td style="text-align:right;">
-2.3417960
</td>
<td style="text-align:right;">
-1.9038976
</td>
</tr>
<tr>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.1282690
</td>
<td style="text-align:right;">
1.6282690
</td>
<td style="text-align:right;">
0.9214369
</td>
<td style="text-align:right;">
0.1985388
</td>
<td style="text-align:right;">
5.8853582
</td>
<td style="text-align:right;">
5.6776743
</td>
<td style="text-align:right;">
-0.2076839
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.8980677
</td>
<td style="text-align:right;">
0.3980677
</td>
<td style="text-align:right;">
-1.1049424
</td>
<td style="text-align:right;">
-0.2965051
</td>
<td style="text-align:right;">
4.7532763
</td>
<td style="text-align:right;">
4.1622951
</td>
<td style="text-align:right;">
-0.5909812
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
-1.3723133
</td>
<td style="text-align:right;">
-1.8723133
</td>
<td style="text-align:right;">
0.9265418
</td>
<td style="text-align:right;">
-0.2643861
</td>
<td style="text-align:right;">
-0.7917697
</td>
<td style="text-align:right;">
-3.5744031
</td>
<td style="text-align:right;">
-2.7826334
</td>
</tr>
</tbody>
</table>

``` r

ptm <- proc.time()
sim_results <- EffectXshift(
  w = w,
  a = a,
  y = y,
  deltas = deltas,
  n_folds = 5,
  num_cores = 6,
  outcome_type = "continuous",
  seed = seed,
  top_n = 1,
  density_classification = TRUE
)
proc.time() - ptm
#>    user  system elapsed 
#>  27.683   1.113 578.345

## marginal effects
k_fold_results <- sim_results$`Effect Modification K-Fold Results`
pooled_results_v <- sim_results$`Effect Modification Region V Pooled Results`
pooled_results_vc <- sim_results$`Effect Modification Region V^c Pooled Results`

## diagnostics to review before interpreting pooled results
selection_diagnostics <- diagnose_selection(sim_results)
positivity_diagnostics <- diagnose_positivity(sim_results)
```

``` r
k_fold_df <- k_fold_results
k_fold_df %>%
  kbl(caption = "K Fold Effect Modification Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
K Fold Effect Modification Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Exposure
</th>
<th style="text-align:right;">
Effect
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower.CI
</th>
<th style="text-align:right;">
Upper.CI
</th>
<th style="text-align:left;">
Modifier
</th>
<th style="text-align:right;">
Fold
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.7579587
</td>
<td style="text-align:right;">
0.1690658
</td>
<td style="text-align:right;">
-1.0893
</td>
<td style="text-align:right;">
-0.4266
</td>
<td style="text-align:left;">
W3 == 0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-1.3879808
</td>
<td style="text-align:right;">
0.3151326
</td>
<td style="text-align:right;">
-2.0056
</td>
<td style="text-align:right;">
-0.7703
</td>
<td style="text-align:left;">
W3 == 1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.5919865
</td>
<td style="text-align:right;">
0.1293933
</td>
<td style="text-align:right;">
-0.8456
</td>
<td style="text-align:right;">
-0.3384
</td>
<td style="text-align:left;">
W3 == 0
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-1.0801027
</td>
<td style="text-align:right;">
0.2593741
</td>
<td style="text-align:right;">
-1.5885
</td>
<td style="text-align:right;">
-0.5717
</td>
<td style="text-align:left;">
W3 == 1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.4721210
</td>
<td style="text-align:right;">
0.1139173
</td>
<td style="text-align:right;">
-0.6954
</td>
<td style="text-align:right;">
-0.2488
</td>
<td style="text-align:left;">
W3 == 0
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-1.1524121
</td>
<td style="text-align:right;">
0.2181057
</td>
<td style="text-align:right;">
-1.5799
</td>
<td style="text-align:right;">
-0.7249
</td>
<td style="text-align:left;">
W3 == 1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.4892525
</td>
<td style="text-align:right;">
0.0949592
</td>
<td style="text-align:right;">
-0.6754
</td>
<td style="text-align:right;">
-0.3031
</td>
<td style="text-align:left;">
W3 == 0
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-1.1872886
</td>
<td style="text-align:right;">
0.2231328
</td>
<td style="text-align:right;">
-1.6246
</td>
<td style="text-align:right;">
-0.7500
</td>
<td style="text-align:left;">
W3 == 1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.2172513
</td>
<td style="text-align:right;">
0.0988232
</td>
<td style="text-align:right;">
-0.4109
</td>
<td style="text-align:right;">
-0.0236
</td>
<td style="text-align:left;">
W3 == 0
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:right;">
-0.7606050
</td>
<td style="text-align:right;">
0.2169844
</td>
<td style="text-align:right;">
-1.1859
</td>
<td style="text-align:right;">
-0.3353
</td>
<td style="text-align:left;">
W3 == 1
</td>
<td style="text-align:right;">
5
</td>
</tr>
</tbody>
</table>

Above the k-fold specific results show that we find the correct
exposure-covariate combination with max effect modification compared to
ground truth. Our estimates are very close to ground-truth as well for a
-0.5 shift in A1 in each level of sex. Above Psi shows the expected
outcome under a shift of -0.5 compared to the observed average outcome.
Variance is the variance derived from the influence function for
stochastic shift interventions. Covariate region is the level of the
covariate that was found. So here, we are looking for the top effect
modification (rank 1), which was found to have exposure A1 in all folds
with effect modification by sex. So in this simple example of 3
exposures and 2 covariates we identify the correct exposure-covariate
combination that has modification and generate estimates that are close
to the truth.

------------------------------------------------------------------------

## For Trialists: Single Binary Treatment

For a randomized trial with one binary treatment, set `rct = TRUE` and
pass one 0/1 treatment column as `a`. The trial mode asks:

> Which baseline covariate-defined subgroup has the largest treatment
> effect difference compared with its complement?

In `rct_type = "ate"` mode, EffectXshift estimates the subject-level
treatment effect *Q*(1, *W*) − *Q*(0, *W*), searches for the region `V`
with the largest differential effect, and reports held-out estimates for
`V`, `V^c`, and the `V - V^c` contrast. If the allocation probability is
known, pass it through `alpha`; if omitted, it is estimated within each
training fold from the observed treatment proportion. For unequal
randomization, the design value is preferred.

Use `target = "risk"` only for a prognostic high-risk subgroup analysis.
That mode finds baseline covariate regions with high held-out outcome
risk; it is not a treatment-effect-modification estimand.

Trialist reporting checklist:

- treatment coding and allocation probability `alpha`
- endpoint and baseline covariates eligible for subgroup discovery
- `n_folds`, `min_obs`, `max_depth`, and `pval_thresh`
- fold-level selected rules, not only the pooled result
- `V`, `V^c`, and `V - V^c` estimates with confidence intervals
- selected-region arm counts and observed outcome summaries from
  `Trial Region Diagnostics`
- whether the analysis was prespecified or exploratory

Current scope: one binary treatment and a marginal treatment
probability. For cluster-randomized, adaptive, crossover, or strongly
stratified designs, the randomization mechanism may need design-specific
handling before trial-grade confirmatory use.

### Fixed-Time Endpoints and Informative Censoring

For a trial endpoint defined at a particular time point *τ*, write the
estimand before running subgroup discovery. For a binary event outcome
this is typically a risk difference,

*E*{1(*T*<sup>1</sup> ≤ *τ*) − 1(*T*<sup>0</sup> ≤ *τ*)},

or the corresponding survival difference. The current `rct = TRUE`
workflow expects that the supplied `y` is already the scalar endpoint to
be analyzed at *τ*. It is appropriate when the endpoint is fully
observed, or when censoring/missingness has already been handled by an
estimator aligned with the trial estimand.

Informative censoring is not handled by coding censored-before-*τ*
subjects as event-free. That instead changes the endpoint and can bias
the fixed-time ATE unless it is the planned composite or
treatment-policy estimand. For a censored time-to-event endpoint, a
future extension should accept event time, event indicator, censoring
time/status, and *τ*; estimate the censoring survival
*G*<sub>*C*</sub>(*t* \| *A*, *W*); form an IPCW/AIPW/TMLE
pseudo-outcome for the fixed-time risk or survival contrast; and then
run the same held-out subgroup discovery on that censoring-adjusted
effect. This is consistent with the ICH E9(R1) emphasis on aligning the
estimand, estimator, and sensitivity analyses (International Council for
Harmonisation of Technical Requirements for Pharmaceuticals for Human
Use 2019) and with targeted-learning approaches for right-censored trial
endpoints (Moore and van der Laan 2009; Brooks et al. 2013).

For a trial report or SAP, also state the intercurrent-event strategy,
whether censoring is assumed conditionally independent given measured
variables, the censoring model/learner, truncation of censoring weights,
and sensitivity analyses for departures from the censoring assumption.

``` r
# A randomized binary treatment whose effect is modified by W1
n <- 1000
W1 <- rbinom(n, 1, 0.5)
W2 <- rnorm(n)
W3 <- rbinom(n, 1, 0.4)
A  <- rbinom(n, 1, 0.5)                       # randomized, alpha = 0.5
tau <- ifelse(W1 == 1, 2, 0)                  # effect modified by W1
Y  <- 1 + 0.5 * W2 - 0.3 * W3 + tau * A + rnorm(n)

rct_results <- EffectXshift(
  w = data.frame(W1 = W1, W2 = W2, W3 = W3),
  a = data.frame(A = A),
  deltas = 0.1,
  y = Y,
  rct = TRUE,
  rct_type = "ate",   # subject-level ATE; use "incps" for an incremental shift
  alpha = 0.5,        # known allocation probability; important if unequal
  n_folds = 2,
  max_depth = 1,
  seed = 4291531
)

# Per-fold discovered regions and effects
rct_results$`Effect Modification K-Fold Results`

# Pooled region V, region V^c, and the V - V^c oracle contrast (with 95% CIs)
rct_results$`Pooled Region Effects`

# Descriptive arm counts and observed outcomes in the selected regions
rct_results$`Trial Region Diagnostics`
```

The `Pooled Region Effects` table reports the effect in `V`, in `V^c`,
and the `V - V^c` contrast, each with a standard error, a Wald 95%
confidence interval, and a p-value. For the data-generating process above
the method recovers the region `{W1 == 1}` with an effect near 2, an
effect near 0 in the complement, and a positive, significant contrast.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/EffectXshift/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/EffectXshift/blob/main/contributing.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/EffectXshift/blob/main/contributing.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `EffectXshift` R package, please cite the following:

------------------------------------------------------------------------

## Related

- [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
  providing an independent implementation of the same core routines for
  the TML estimation procedure and statistical methodology as is made
  available here, through reliance on a unified interface for Targeted
  Learning provided by the [`tmle3`](https://github.com/tlverse/tmle3)
  engine of the [`tlverse` ecosystem](https://github.com/tlverse).

- [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
  providing facilities to estimate the causal effect of stochastic
  treatment regimes in the mediation setting, including classical (IPW)
  and augmented double robust (one-step) estimators. This is an
  implementation of the methodology explored by Dı́az and Hejazi (2020).

- [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
  package for estimating the conditional density treatment mechanism
  component of this parameter based on using the [highly adaptive
  lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, Phillips,
  et al. 2022; Hejazi, Coyle, and van der Laan 2020) in combination with
  a pooled hazard regression. This package implements a variant of the
  approach advocated by Dı́az and van der Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through NIH grant
P42ES004705 from NIEHS

------------------------------------------------------------------------

## License

© 2020-2022 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2020-2024 David B. McCoy
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-athey2016recursive" class="csl-entry">

Athey, Susan, and Guido Imbens. 2016. “Recursive Partitioning for
Heterogeneous Causal Effects.” *Proceedings of the National Academy of
Sciences* 113 (27): 7353–60.
<https://doi.org/10.1073/pnas.1510489113>.

</div>

<div id="ref-brooks2013targeted" class="csl-entry">

Brooks, Jordan C, Mark J van der Laan, Daniel E Singer, and Alan S Go.
2013. “Targeted Minimum Loss-Based Estimation of Causal Effects in
Right-Censored Survival Data with Time-Dependent Covariates: Warfarin,
Stroke, and Death in Atrial Fibrillation.” *Journal of Causal Inference*
1 (2): 235–54. <https://doi.org/10.1515/jci-2013-0001>.

</div>

<div id="ref-coyle-sl3-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. “<span class="nocase">sl3</span>: Modern Machine
Learning Pipelines for Super Learning.”
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle-hal9001-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. “<span class="nocase">hal9001</span>: The
Scalable Highly Adaptive Lasso.”
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal" class="csl-entry">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population" class="csl-entry">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-haneuse2013estimation" class="csl-entry">

Haneuse, Sebastian, and Andrea Rotnitzky. 2013. “Estimation of the
Effect of Interventions That Modify the Received Treatment.” *Statistics
in Medicine* 32 (30): 5260–77.

</div>

<div id="ref-ich2019e9r1" class="csl-entry">

International Council for Harmonisation of Technical Requirements for
Pharmaceuticals for Human Use. 2019. “ICH E9(R1) Addendum on Estimands
and Sensitivity Analysis in Clinical Trials to the Guideline on
Statistical Principles for Clinical Trials.”
<https://database.ich.org/sites/default/files/E9-R1_Step4_Guideline_2019_1203.pdf>.

</div>

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020.
“<span class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso
Regression in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-keil2020quantile" class="csl-entry">

Keil, Alexander P, Jessie P Buckley, Katie M O'Brien, Kelly K Ferguson,
Shanshan Zhao, and Alexandra J White. 2020. “A Quantile-Based
G-Computation Approach to Addressing the Effects of Exposure Mixtures.”
*Environmental Health Perspectives* 128 (4): 047004.
<https://doi.org/10.1289/EHP5838>.

</div>

<div id="ref-kunzel2019metalearners" class="csl-entry">

Künzel, Sören R, Jasjeet S Sekhon, Peter J Bickel, and Bin Yu. 2019.
“Metalearners for Estimating Heterogeneous Treatment Effects Using
Machine Learning.” *Proceedings of the National Academy of Sciences* 116
(10): 4156–65. <https://doi.org/10.1073/pnas.1804597116>.

</div>

<div id="ref-moore2009increasing" class="csl-entry">

Moore, Kelly L, and Mark J van der Laan. 2009. “Increasing Power in
Randomized Trials with Right Censored Outcomes Through Covariate
Adjustment.” *Journal of Biopharmaceutical Statistics* 19 (6): 1099–1131.
<https://doi.org/10.1080/10543400903243017>.

</div>

<div id="ref-wager2018estimation" class="csl-entry">

Wager, Stefan, and Susan Athey. 2018. “Estimation and Inference of
Heterogeneous Treatment Effects Using Random Forests.” *Journal of the
American Statistical Association* 113 (523): 1228–42.
<https://doi.org/10.1080/01621459.2017.1319839>.

</div>

</div>
