# ---------------------------------------------------------------------------
# Diagnostic: is the mixture path's pooled undercoverage driven by BIAS
# (finite-sample, shrinks with n) or by VARIANCE UNDERESTIMATION (reported SE
# too small -- a calibration bug independent of n)?
#
# For each replicate we record the pooled point estimate AND its reported SE for
# region V, region V^c, and the V - V^c contrast. We then compare the *mean
# reported SE* to the *empirical SD of the estimates across replicates*:
#   ratio = mean(reported SE) / sd(estimates)
#   ratio ~ 1   -> variance estimator is calibrated; undercoverage is bias.
#   ratio << 1  -> reported SEs are too small; a variance bug.
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

set.seed(20240611)

n_reps  <- as.integer(Sys.getenv("DIAG_REPS", "30"))
n_obs   <- as.integer(Sys.getenv("DIAG_N", "1000"))
n_folds <- as.integer(Sys.getenv("DIAG_FOLDS", "2"))
delta   <- 1
true_v  <- 3 * delta
true_vc <- 1 * delta
true_con <- true_v - true_vc

sim_one <- function() {
  W1 <- rbinom(n_obs, 1, 0.5); W2 <- rnorm(n_obs); W3 <- rnorm(n_obs)
  A1 <- rnorm(n_obs, 0.3 * W2); A2 <- rnorm(n_obs, 0.2 * W3); A3 <- rnorm(n_obs)
  Y  <- 1 + A1 + 0.5 * A2 + 0.2 * A3 + 2 * A1 * W1 + 0.4 * W2 + rnorm(n_obs)
  w <- data.frame(W1 = W1, W2 = W2, W3 = W3)
  a <- data.frame(A1 = A1, A2 = A2, A3 = A3)
  ml <- list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new())

  dens_class <- as.logical(Sys.getenv("DIAG_DENSITY_CLASS", "TRUE"))
  # With direct density estimation use the package's default density Super
  # Learner (g_learner = NULL builds it), since a regression stack is not a
  # density estimator.
  g_lrn <- if (dens_class) ml else NULL

  res <- tryCatch(EffectXshift(
    w = w, a = a, y = Y, deltas = list(A1 = delta, A2 = delta, A3 = delta),
    mu_learner = ml, g_learner = g_lrn, n_folds = n_folds, parallel = FALSE, seed = 1,
    rct = FALSE, density_classification = dens_class, min_obs = 40, max_depth = 1, top_n = 1
  ), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  v  <- res[["Effect Modification Region V Pooled Results"]]
  vc <- res[["Effect Modification Region V^c Pooled Results"]]
  if (is.null(v) || is.null(vc)) return(NULL)
  data.frame(
    v_est = v$Psi, v_se = v$SE,
    vc_est = vc$Psi, vc_se = vc$SE,
    con_est = v$Psi - vc$Psi, con_se = sqrt(v$SE^2 + vc$SE^2)
  )
}

r <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
  if (i %% 5 == 0) message("rep ", i, "/", n_reps); sim_one()
}))

report <- function(name, est, se, truth) {
  emp_sd <- sd(est); mean_se <- mean(se)
  cov <- mean(est - 1.96 * se <= truth & est + 1.96 * se >= truth)
  cat(sprintf("%-10s  est=%.3f (truth %.1f)  empSD=%.3f  meanSE=%.3f  SE/SD=%.2f  cov=%.2f\n",
              name, mean(est), truth, emp_sd, mean_se, mean_se / emp_sd, cov))
}

cat("\n===== Mixture variance calibration (", nrow(r), "reps, n =", n_obs,
    ", folds =", n_folds, ") =====\n")
report("Region V",  r$v_est,  r$v_se,  true_v)
report("Region V^c", r$vc_est, r$vc_se, true_vc)
report("Contrast",  r$con_est, r$con_se, true_con)
cat("\nSE/SD ~ 1 => calibrated (undercoverage is bias). SE/SD << 1 => SEs too small.\n")

# ---------------------------------------------------------------------------
# Reference results (ranger(100) + glm, 2 folds; "Contrast" = V - V^c, truth 2).
#
#                       bias    SE/SD   coverage
#   PRE-FIX  n=1000    -0.41    0.16      0.20      (no-shift density ratio == 1)
#   POST-FIX n=1000    -0.21    0.48      0.63
#   POST-FIX n=2500    -0.00    0.38      0.50
#
# Reading the results:
#  * The density-ratio fix (no-shift clever covariate = classifier odds, not 1)
#    restores debiasing: the bias now goes to ~0 as n grows, and the egregious
#    variance underestimation (SE/SD 0.16) roughly triples toward calibration.
#  * A residual ~2x SE underestimation remains and does NOT vanish with n in
#    this regime: the empirical SD shrinks slower than the 1/sqrt(n) rate of the
#    influence-function SE, because flexible-ML nuisance error (the ranger
#    outcome/density fits, which converge at a nonparametric rate) is not
#    captured by the first-order EIF with only 2 cross-fitting folds.
#  * Practical guidance: for trustworthy CIs use the paper's regime -- a richer
#    Super Learner library, >= 5 folds, and larger n -- which speeds nuisance
#    convergence so the first-order variance becomes calibrated. The RCT path
#    (known propensity, AIPW) is unaffected and stays at ~95% coverage.
#
# Methodological-lever test (does increasing folds alone fix calibration?):
#                       contrast SE/SD   contrast coverage
#   2 folds (n=1000)        0.48              0.63
#   5 folds (n=1000)        0.56              0.60      (15 reps; noisy)
#   -> More folds alone gives only a marginal change. The bottleneck is the
#      learner's bias / slow convergence, not the cross-fitting count, so the
#      real remedy is a faster-converging nuisance library (e.g. HAL with
#      undersmoothing) or a resampling-based standard error. The default n_folds
#      was raised to 5 as better cross-fitting practice (matching the paper), but
#      it is NOT on its own a cure for the finite-sample CI undercoverage.
# ---------------------------------------------------------------------------
