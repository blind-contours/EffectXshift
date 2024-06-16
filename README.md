
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`EffectXshift`

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/EffectXshift/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/EffectXShift/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/EffectXshift/master.svg)](https://codecov.io/github/blind-contours/EffectXShift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/EffectXshift)](https://www.r-pkg.org/pkg/EffectXshift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/EffectXShift)](https://CRAN.R-project.org/package=EffectXshift)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/EffectXShift)](https://CRAN.R-project.org/package=EffectXshift)
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
an outcome. We define effect modification as the counterfactual mean of
the outcome under stochastic interventions of an exposure in the mixture
in a subregion of the covariate space compared to the complementary
region of that space. Stochastic interventions or exposure changes
depend on naturally observed values, as described in past literature
(Dı́az and van der Laan 2012; Haneuse and Rotnitzky 2013).

Our target parameter is:

$$
\begin{align*}
&\underset{A_i \in \boldsymbol{A}, \, V \subseteq W}{\text{arg max}} \\
&\quad E\left[\, E\left[Y \,|\, A_i + \delta_i, \, W_j \in V, \, W_{\setminus j}\right] \right. \\
&\quad \left. - E\left[Y \,|\,A_i + \delta_i, \, W_j \in V^c, \, W_{\setminus j}\right] \, \right].
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

The package ensures robustness by employing a k-fold cross-validation
framework. This framework helps in estimating a data-adaptive parameter,
which is the stochastic shift target parameters for the exposure
identified as having maximal effect modification. The process begins by
partitioning the data into parameter-generating and estimation samples.
In the parameter-generating sample, we identify our exposure-covariate
region with maximum effect modification using machine learning
g-computation framework. In the estimation sample we then estimate our
stochastic shift target parameter in levels of the discovered covariate
using the doubly robust estimator TMLE to ensure asymptotic efficiency
which allows us to construct confidence intervals for our estimates
(unlike the g-comp method).

By using EffectXshift, users get access to a tool that offers both
k-fold specific and aggregated results for the maximal effect
modification, ensuring that researchers can glean the most information
from their data. For a more in-depth exploration, there’s an
accompanying vignette.

## Inputs

To utilize the package, users need to provide vectors for exposures,
covariates, and outcomes. They also specify the respective $\delta$ for
each exposure (indicating the degree of shift). The `top_n` parameter
defines the top number of effect modification exposure-covariate pairs.
A detailed guide is provided in the vignette. With these inputs,
`EffectXshift` processes the data and delivers tables showcasing
fold-specific results and aggregated outcomes, allowing users to glean
insights effectively.

`EffectXshift` also incorporates features from the `sl3` package (Coyle,
Hejazi, Malenica, et al. 2022), facilitating ensemble machine learning
in the estimation process. If the user does not specify any stack
parameters, `EffectXshift` will automatically create an ensemble of
machine learning algorithms that strike a balance between flexibility
and computational efficiency.

------------------------------------------------------------------------

## Installation

*Note:* Because the `EffectXshift` package (currently) depends on `sl3`
that allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `EffectXshift` package is not
available on CRAN and must be downloaded here.

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `EffectXshift`

``` r
remotes::install_github("blind-contours/EffectXshift@main")
```

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
n <- 1000 # Number of observations

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
#> 1     0     0     0         -0.523
#> 2     0     0     1         -0.501
#> 3     0     1     0         -0.449
#> 4     0     1     1         -0.893
#> 5     1     0     0         -1.83 
#> 6     1     0     1         -1.58 
#> 7     1     1     0         -1.43 
#> 8     1     1     1         -1.35

data_shifted$Y_diff <- data_shifted$Y_shifted - data_shifted$Y

# Aggregate differences by levels of Sex
avg_diff_W3 <- aggregate(Y_diff ~ W3, data = data_shifted, mean)

# Aggregate differences by levels of W2
avg_diff_W2 <- aggregate(Y_diff ~ W2, data = data_shifted, mean)

# Aggregate differences by levels of W1
avg_diff_W1 <- aggregate(Y_diff ~ W1, data = data_shifted, mean)

print(avg_diff_W3)
#>   W3     Y_diff
#> 1  0 -0.5845702
#> 2  1 -1.5486749
print(avg_diff_W2)
#>   W2    Y_diff
#> 1  0 -1.108985
#> 2  1 -1.051339
print(avg_diff_W1)
#>   W1    Y_diff
#> 1  0 -1.081664
#> 2  1 -1.080522
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
1.8251942
</td>
<td style="text-align:right;">
1.3251942
</td>
<td style="text-align:right;">
-0.0451165
</td>
<td style="text-align:right;">
0.3677981
</td>
<td style="text-align:right;">
9.6444070
</td>
<td style="text-align:right;">
6.4087359
</td>
<td style="text-align:right;">
-3.2356710
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
1.4355733
</td>
<td style="text-align:right;">
0.9355733
</td>
<td style="text-align:right;">
-0.5961851
</td>
<td style="text-align:right;">
-0.1900499
</td>
<td style="text-align:right;">
6.7776261
</td>
<td style="text-align:right;">
3.9907013
</td>
<td style="text-align:right;">
-2.7869247
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
-0.3597802
</td>
<td style="text-align:right;">
-0.8597802
</td>
<td style="text-align:right;">
0.2278567
</td>
<td style="text-align:right;">
0.4210997
</td>
<td style="text-align:right;">
0.8195058
</td>
<td style="text-align:right;">
0.2503849
</td>
<td style="text-align:right;">
-0.5691209
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
0
</td>
<td style="text-align:right;">
-0.1521016
</td>
<td style="text-align:right;">
-0.6521016
</td>
<td style="text-align:right;">
-0.5440319
</td>
<td style="text-align:right;">
-0.0723255
</td>
<td style="text-align:right;">
2.3289856
</td>
<td style="text-align:right;">
1.2345317
</td>
<td style="text-align:right;">
-1.0944539
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
0
</td>
<td style="text-align:right;">
0.5522047
</td>
<td style="text-align:right;">
0.0522047
</td>
<td style="text-align:right;">
3.3883258
</td>
<td style="text-align:right;">
1.1761818
</td>
<td style="text-align:right;">
7.5480306
</td>
<td style="text-align:right;">
2.8384768
</td>
<td style="text-align:right;">
-4.7095538
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
0
</td>
<td style="text-align:right;">
1.2299938
</td>
<td style="text-align:right;">
0.7299938
</td>
<td style="text-align:right;">
2.4612301
</td>
<td style="text-align:right;">
0.5439224
</td>
<td style="text-align:right;">
5.7973901
</td>
<td style="text-align:right;">
5.1828397
</td>
<td style="text-align:right;">
-0.6145505
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
#>   6.770   0.454  87.132

## marginal effects
k_fold_results <- sim_results$`Effect Modification K-Fold Results`
pooled_results_v <- sim_results$`Effect Modification Region V Pooled Results`
pooled_results_vc <- sim_results$`Effect Modification Region V^c Pooled Results`
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
-0.4607403
</td>
<td style="text-align:right;">
0.1716025
</td>
<td style="text-align:right;">
-0.7971
</td>
<td style="text-align:right;">
-0.1244
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
-0.9372819
</td>
<td style="text-align:right;">
0.3293086
</td>
<td style="text-align:right;">
-1.5827
</td>
<td style="text-align:right;">
-0.2918
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
-0.4312649
</td>
<td style="text-align:right;">
0.1520413
</td>
<td style="text-align:right;">
-0.7293
</td>
<td style="text-align:right;">
-0.1333
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
-0.9259123
</td>
<td style="text-align:right;">
0.4038112
</td>
<td style="text-align:right;">
-1.7174
</td>
<td style="text-align:right;">
-0.1345
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
-0.9851225
</td>
<td style="text-align:right;">
0.1663141
</td>
<td style="text-align:right;">
-1.3111
</td>
<td style="text-align:right;">
-0.6592
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
-1.7650811
</td>
<td style="text-align:right;">
0.2966523
</td>
<td style="text-align:right;">
-2.3465
</td>
<td style="text-align:right;">
-1.1837
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
-0.3049075
</td>
<td style="text-align:right;">
0.2336208
</td>
<td style="text-align:right;">
-0.7628
</td>
<td style="text-align:right;">
0.1530
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
-0.8764155
</td>
<td style="text-align:right;">
0.3151982
</td>
<td style="text-align:right;">
-1.4942
</td>
<td style="text-align:right;">
-0.2586
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
-0.1232804
</td>
<td style="text-align:right;">
0.1321311
</td>
<td style="text-align:right;">
-0.3823
</td>
<td style="text-align:right;">
0.1357
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
-0.7886600
</td>
<td style="text-align:right;">
0.5435316
</td>
<td style="text-align:right;">
-1.8540
</td>
<td style="text-align:right;">
0.2766
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

The consistency of our results means we can look at the pooled results
which pools estimates for the top effect modifier in each level.

``` r
pooled_results_df <- rbind(pooled_results_v, pooled_results_vc)
pooled_results_df %>%
  kbl(caption = "Pooled TMLE Results by Region") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Pooled TMLE Results by Region
</caption>
<thead>
<tr>
<th style="text-align:left;">
Condition
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:left;">
Type
</th>
<th style="text-align:left;">
Variables
</th>
<th style="text-align:right;">
N
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
v
</td>
<td style="text-align:right;">
0.1985276
</td>
<td style="text-align:right;">
0.0048664
</td>
<td style="text-align:right;">
0.0697597
</td>
<td style="text-align:right;">
0.0618
</td>
<td style="text-align:right;">
0.3353
</td>
<td style="text-align:right;">
0.0044289
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
v
</td>
<td style="text-align:right;">
485
</td>
</tr>
<tr>
<td style="text-align:left;">
vc
</td>
<td style="text-align:right;">
-1.4408420
</td>
<td style="text-align:right;">
0.0220128
</td>
<td style="text-align:right;">
0.1483672
</td>
<td style="text-align:right;">
-1.7316
</td>
<td style="text-align:right;">
-1.1500
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
vc
</td>
<td style="text-align:right;">
515
</td>
</tr>
</tbody>
</table>

This shows the pooled results for the region v and its compliment.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/EffectXshift/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20EffectXshift/main/contributing.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/EffectXshift/blob/master/CONTRIBUTING.md)
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

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020.
“<span class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso
Regression in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

</div>
