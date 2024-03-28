
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`EffectXshift` <img src="man/figures/man/figures/EffectXshift_logo.png" style="float:right; height:200px;">

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
the outcome under stochastic interventions in a subregion of the
exposure space compared to the complementary region of that space.
Stochastic interventions or exposure changes depend on naturally
observed values, as described in past literature (Dı́az and van der Laan
2012; Haneuse and Rotnitzky 2013).

Our target parameter is:

$$
\begin{align*}
&\underset{A_i \in \boldsymbol{A}, \, V \subseteq W}{\text{arg max}} \\
&\quad E\left[\, E\left[Y \,|\, A_i + \delta_i, \, W_j \in V, \, W_{\setminus j}\right] \right. \\
&\quad \left. - E\left[Y \,|\,A_i + \delta_i, \, W_j \in V^c, \, W_{\setminus j}\right] \, \right].
\end{align*}
$$ 

This finds the exposure-covariate region combination where the effect
of an intervention is maximally diffrent.

`EffectXshift` builds on work described in (McCoy et al. 2023). However
instead of identifying interactions through an semi-parametric
definition of an F-statistics and then estimating our interaction target
parameter using CV-TMLE pooled across exposure sets, we provide a more
streamlined approach. In this package, we identify interactions through
g-computation first - evaluating the expected outcome under joint shift
compared to the sum of individual shifts using Super Learner. We then
rank these estimates as the highest sets of synergy and antagonism. We
then use CV-TMLE and pool within the ranks.

The package ensures robustness by employing a k-fold cross-validation
framework. This framework helps in estimating a data-adaptive parameter,
which is the stochastic shift target parameters for the exposure sets
identified as having synergy or antagonism. The process begins by
partitioning the data into parameter-generating and estimation samples.
In the parameter-generating sample, we identify our ranks of
antagonistic and synergistic exposure sets through a machine learning
g-computation framework. In the estimation sample we then estimate our
interaction target parameter using the doubly robust estimator TMLE to
ensure asymptotic efficiency which allows us to construct confidence
intervals for our estimates (unlike the g-comp method).

By using EffectXshift, users get access to a tool that offers both
k-fold specific and aggregated results for the top synergistic and
antagonistic relationships, ensuring that researchers can glean the most
information from their data. For a more in-depth exploration, there’s an
accompanying vignette.

To utilize the package, users need to provide vectors for exposures,
covariates, and outcomes. They also specify the respective $\delta$ for
each exposure (indicating the degree of shift) and if this delta should
be adaptive in response to positivity violations. The `top_n` parameter
defines the top number of synergistic, antagonistic, positive and
negative ranked impacts to estiamte. A detailed guide is provided in the
vignette. With these inputs, `EffectXshift` processes the data and
delivers tables showcasing fold-specific results and aggregated
outcomes, allowing users to glean insights effectively.

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

There are many depedencies for `EffectXshift` so it’s easier to break up
installation of the various packages to ensure proper installation.

First install the basis estimators used in the data-adaptive variable
discovery of the exposure and covariate space:

``` r
install.packages("earth")
install.packages("hal9001")
```

`EffectXshift` uses the `sl3` package to build ensemble machine learners
for each nuisance parameter. We have to install off the development
branch, first download these two packages for `sl3`

``` r
install.packages(c("ranger", "arm", "xgboost", "nnls"))
```

Now install `sl3` on devel:

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `EffectXshift`

``` r
remotes::install_github("blind-contours/EffectXshift@main")
```

`EffectXshift` has some other miscellaneous dependencies that are used
in the examples as well as in the plotting functions.

``` r
install.packages(c("kableExtra", "hrbrthemes", "viridis"))
```

------------------------------------------------------------------------

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
set.seed(123)  # Ensure reproducibility

n <- 1000  # Number of observations

# Generate binary covariates
Sex <- rbinom(n, 1, 0.5)  # 0 for one sex, 1 for another
W2 <- rbinom(n, 1, 0.5)  # Additional binary covariate

# Generate exposures
A1 <- rnorm(n)  # Continuous exposure with significant interaction with Sex
A2 <- rnorm(n)  # Continuous exposure without significant interaction
A3 <- rnorm(n)  # Continuous exposure without significant interaction

# Define effect sizes
beta_A1 <- 1  # Base effect of A1 on Y
beta_A2 <- 0.5  # Effect of A2 on Y
beta_A3 <- 0.2  # Effect of A3 on Y
beta_Sex <- 0.5  # Effect of Sex on Y
beta_W2 <- -0.3  # Effect of W2 on Y
interaction_A1_Sex <- 2  # Interaction effect between A1 and Sex, ensuring maximal impact difference

# Simulate outcome Y
Y <- 2 + beta_A1*A1 + beta_A2*A2 + beta_A3*A3 + beta_Sex*Sex + beta_W2*W2 +
     interaction_A1_Sex*A1*Sex + rnorm(n)  # Include noise


# Assume an intervention that reduces A1 by a fixed amount for the shift analysis
shift_amount <- -0.5
A1_shifted <- A1 + shift_amount

# Recalculate Y assuming the shift in A1
Y_shifted <- 2 + beta_A1*A1_shifted + beta_A2*A2 + beta_A3*A3 + beta_Sex*Sex + beta_W2*W2 +
             interaction_A1_Sex*A1_shifted*Sex + rnorm(n)

# Data frame with shifted outcomes
data_shifted <- data.frame(Sex, W2, A1, A1_shifted, A2, A3, Y, Y_shifted)

effect_summary <- data_shifted %>%
  mutate(Difference = Y_shifted - Y) %>%
  group_by(Sex, W2) %>%
  summarise(
    Avg_Difference = mean(Difference),
    .groups = 'drop'
  )

effect_summary
#> # A tibble: 4 × 3
#>     Sex    W2 Avg_Difference
#>   <int> <int>          <dbl>
#> 1     0     0         -0.466
#> 2     0     1         -0.367
#> 3     1     0         -1.40 
#> 4     1     1         -1.52

data_shifted$Y_diff = data_shifted$Y_shifted - data_shifted$Y

# Aggregate differences by levels of Sex
avg_diff_sex <- aggregate(Y_diff ~ Sex, data = data_shifted, mean)

# Aggregate differences by levels of W2
avg_diff_W2 <- aggregate(Y_diff ~ W2, data = data_shifted, mean)

avg_diff_sex
#>   Sex     Y_diff
#> 1   0 -0.4138098
#> 2   1 -1.4534507
avg_diff_W2
#>   W2     Y_diff
#> 1  0 -0.9507121
#> 2  1 -0.9016993
```

``` r
w <- data_shifted[, c("W2", "Sex")]
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
Sex
</th>
<th style="text-align:right;">
W2
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
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
-0.9957987
</td>
<td style="text-align:right;">
-1.4957987
</td>
<td style="text-align:right;">
-0.5116037
</td>
<td style="text-align:right;">
-0.1503075
</td>
<td style="text-align:right;">
0.9148877
</td>
<td style="text-align:right;">
-0.2758360
</td>
<td style="text-align:right;">
-1.190724
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
-1.0399550
</td>
<td style="text-align:right;">
-1.5399550
</td>
<td style="text-align:right;">
0.2369379
</td>
<td style="text-align:right;">
-0.3277571
</td>
<td style="text-align:right;">
-0.2168344
</td>
<td style="text-align:right;">
-1.2393542
</td>
<td style="text-align:right;">
-1.022520
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
-0.0179802
</td>
<td style="text-align:right;">
-0.5179802
</td>
<td style="text-align:right;">
-0.5415892
</td>
<td style="text-align:right;">
-1.4481653
</td>
<td style="text-align:right;">
2.0925963
</td>
<td style="text-align:right;">
-0.2253574
</td>
<td style="text-align:right;">
-2.317954
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
-0.1321751
</td>
<td style="text-align:right;">
-0.6321751
</td>
<td style="text-align:right;">
1.2192276
</td>
<td style="text-align:right;">
-0.6972846
</td>
<td style="text-align:right;">
0.9894737
</td>
<td style="text-align:right;">
2.2546501
</td>
<td style="text-align:right;">
1.265176
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
-2.5493428
</td>
<td style="text-align:right;">
-3.0493428
</td>
<td style="text-align:right;">
0.1741359
</td>
<td style="text-align:right;">
2.5984902
</td>
<td style="text-align:right;">
-6.8673719
</td>
<td style="text-align:right;">
-5.4250711
</td>
<td style="text-align:right;">
1.442301
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
1.0405735
</td>
<td style="text-align:right;">
0.5405735
</td>
<td style="text-align:right;">
-0.6152683
</td>
<td style="text-align:right;">
-0.0374150
</td>
<td style="text-align:right;">
4.9307824
</td>
<td style="text-align:right;">
2.5605873
</td>
<td style="text-align:right;">
-2.370195
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
  delta = deltas,
  n_folds = 3,
  num_cores = 6,
  outcome_type = "continuous",
  seed = seed,
  top_n = 1
)
#> 
#> Iter: 1 fn: 479.5067  Pars:  0.999994748 0.000005252
#> Iter: 2 fn: 479.5067  Pars:  0.99999858 0.00000142
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 482.2148  Pars:  0.99998925 0.00001075
#> Iter: 2 fn: 482.2148  Pars:  0.999997737 0.000002263
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 491.5928  Pars:  0.00001227 0.99998773
#> Iter: 2 fn: 491.5928  Pars:  0.000005372 0.999994628
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 468.9115  Pars:  0.55287 0.44713
#> Iter: 2 fn: 468.9106  Pars:  0.50547 0.49453
#> Iter: 3 fn: 468.9106  Pars:  0.50547 0.49453
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 477.9934  Pars:  0.999993554 0.000006445
#> Iter: 2 fn: 477.9934  Pars:  0.999998444 0.000001556
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 485.8838  Pars:  0.99997075 0.00002925
#> Iter: 2 fn: 485.8838  Pars:  0.999991334 0.000008666
#> solnp--> Completed in 2 iterations
proc.time() - ptm
#>    user  system elapsed 
#>  12.117   0.727 183.695

## marginal effects
k_fold_results <- sim_results$`Effect Modification K-Fold Results`
pooled_results <- sim_results$`Effect Modification Pooled Results`
```

``` r
k_fold_df <- do.call(rbind, k_fold_results)
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
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
<th style="text-align:left;">
Covariate_region
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fold : 1 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 1
</td>
<td style="text-align:right;">
-0.5917534
</td>
<td style="text-align:right;">
0.0098075
</td>
<td style="text-align:right;">
0.0990326
</td>
<td style="text-align:right;">
-0.7859
</td>
<td style="text-align:right;">
-0.3977
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
170
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fold : 1 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 2
</td>
<td style="text-align:right;">
-1.4752810
</td>
<td style="text-align:right;">
0.0555699
</td>
<td style="text-align:right;">
0.2357326
</td>
<td style="text-align:right;">
-1.9373
</td>
<td style="text-align:right;">
-1.0133
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
164
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 1
</td>
</tr>
<tr>
<td style="text-align:left;">
Fold : 2 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 1
</td>
<td style="text-align:right;">
-0.5428111
</td>
<td style="text-align:right;">
0.0094160
</td>
<td style="text-align:right;">
0.0970359
</td>
<td style="text-align:right;">
-0.7330
</td>
<td style="text-align:right;">
-0.3526
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fold : 2 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 2
</td>
<td style="text-align:right;">
-1.5179830
</td>
<td style="text-align:right;">
0.0602319
</td>
<td style="text-align:right;">
0.2454219
</td>
<td style="text-align:right;">
-1.9990
</td>
<td style="text-align:right;">
-1.0370
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
172
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 1
</td>
</tr>
<tr>
<td style="text-align:left;">
Fold : 3 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 1
</td>
<td style="text-align:right;">
-0.5477520
</td>
<td style="text-align:right;">
0.0114864
</td>
<td style="text-align:right;">
0.1071745
</td>
<td style="text-align:right;">
-0.7578
</td>
<td style="text-align:right;">
-0.3377
</td>
<td style="text-align:right;">
3e-07
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
176
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fold : 3 \| Rank : 1 \| Exposure : A1 \| Modifier : Sex \| Level : 2
</td>
<td style="text-align:right;">
-1.4742615
</td>
<td style="text-align:right;">
0.0704286
</td>
<td style="text-align:right;">
0.2653839
</td>
<td style="text-align:right;">
-1.9944
</td>
<td style="text-align:right;">
-0.9541
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
157
</td>
<td style="text-align:right;">
-0.5
</td>
<td style="text-align:left;">
Sex == 1
</td>
</tr>
</tbody>
</table>

Above the k-fold specific results show that we find the correct
exposure-covariate combination with max effect modification compared to
ground truth. Our estimates are very close to ground-truth as well for a
-0.5 shift in A1 in each level of sex.

The consistency of our results means we can look at the pooled results
which pools estimates for the top effect modifier in each level.

``` r
pooled_results_df <- do.call(rbind, pooled_results)
pooled_results_df %>%
  kbl(caption = "Pooled TMLE Results by Level") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Pooled TMLE Results by Level
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
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:right;">
-0.5632749
</td>
<td style="text-align:right;">
0.0038986
</td>
<td style="text-align:right;">
0.0624389
</td>
<td style="text-align:right;">
-0.6857
</td>
<td style="text-align:right;">
-0.4409
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
507
</td>
<td style="text-align:right;">
-0.5
</td>
</tr>
<tr>
<td style="text-align:left;">
1_2
</td>
<td style="text-align:right;">
-1.5033988
</td>
<td style="text-align:right;">
0.0211669
</td>
<td style="text-align:right;">
0.1454886
</td>
<td style="text-align:right;">
-1.7886
</td>
<td style="text-align:right;">
-1.2182
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
493
</td>
<td style="text-align:right;">
-0.5
</td>
</tr>
</tbody>
</table>

This shows that in rank 1 level 1 (1_1) the pooled result which pools
our fold estimates for level 1 for the top modifier. 1_2 is the rank 1
and second level of modifier. These correspond to sex = 0, and sex = 1
in this simulated case.

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

The development of this software was supported in part through grants
from the

------------------------------------------------------------------------

## License

© 2020-2022 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2020-2022 David B. McCoy
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

<div id="ref-mccoy2023semiparametric" class="csl-entry">

McCoy, David B., Alan E. Hubbard, Alejandro Schuler, and Mark J. van der
Laan. 2023. “Semi-Parametric Identification and Estimation of
Interaction and Effect Modification in Mixed Exposures Using Stochastic
Interventions.” <https://arxiv.org/abs/2305.01849>.

</div>

</div>
