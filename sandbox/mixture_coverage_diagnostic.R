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

n_reps <- as.integer(Sys.getenv("DIAG_REPS", "30"))
n_obs  <- as.integer(Sys.getenv("DIAG_N", "1000"))
delta  <- 1
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
    mu_learner = ml, g_learner = g_lrn, n_folds = 2, parallel = FALSE, seed = 1,
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

cat("\n===== Mixture variance calibration (", nrow(r), "reps, n =", n_obs, ") =====\n")
report("Region V",  r$v_est,  r$v_se,  true_v)
report("Region V^c", r$vc_est, r$vc_se, true_vc)
report("Contrast",  r$con_est, r$con_se, true_con)
cat("\nSE/SD ~ 1 => calibrated (undercoverage is bias). SE/SD << 1 => SEs too small.\n")
