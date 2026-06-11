# ---------------------------------------------------------------------------
# Simulation validation for the RCT workflow (EffectXshift(rct = TRUE)).
#
# DGP: randomized binary treatment A (alpha = 0.5) whose effect is modified by a
# binary covariate W1. The conditional ATE is 2 when W1 == 1 and 0 when W1 == 0.
# Therefore the oracle region V = {W1 == 1} has true effect 2, V^c has true
# effect 0, and the oracle V - V^c contrast equals 2.
#
# We report, over many replicates:
#   * detection rate  : how often the discovered region rule references W1
#   * bias / SE / RMSE : for region V, region V^c, and the contrast
#   * CI coverage      : Wald 95% CI coverage of each true value
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

set.seed(20240611)

n_reps   <- 100
n_obs    <- 1000
true_v   <- 2 # true ATE in region V = {W1 == 1}
true_vc  <- 0 # true ATE in region V^c
true_con <- true_v - true_vc # true oracle contrast

sim_one <- function() {
  W1 <- rbinom(n_obs, 1, 0.5)
  W2 <- rnorm(n_obs)
  W3 <- rbinom(n_obs, 1, 0.4)
  A  <- rbinom(n_obs, 1, 0.5)
  tau <- ifelse(W1 == 1, 2, 0)
  Y  <- 1 + 0.5 * W2 - 0.3 * W3 + tau * A + rnorm(n_obs)

  w <- data.frame(W1 = W1, W2 = W2, W3 = W3)
  a <- data.frame(A = A)
  mu_learner <- list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new())

  res <- tryCatch(
    EffectXshift(
      w = w, a = a, y = Y, deltas = 0.1,
      mu_learner = mu_learner, g_learner = mu_learner,
      n_folds = 2, parallel = FALSE, seed = 1,
      rct = TRUE, rct_type = "ate",
      min_obs = 40, max_depth = 1, pval_thresh = 0.2
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  pooled <- res[["Pooled Region Effects"]]
  kfold  <- res[["Effect Modification K-Fold Results"]]
  get <- function(region) pooled[pooled$Region == region, ]
  v  <- get("V"); vc <- get("V^c"); con <- get("V - V^c")

  data.frame(
    detected = any(grepl("W1", kfold$Rule)),
    v_est = v$Psi,  v_cov = (v$`Lower CI`  <= true_v   & v$`Upper CI`  >= true_v),
    vc_est = vc$Psi, vc_cov = (vc$`Lower CI` <= true_vc  & vc$`Upper CI` >= true_vc),
    con_est = con$Psi, con_se = con$SE,
    con_cov = (con$`Lower CI` <= true_con & con$`Upper CI` >= true_con)
  )
}

results <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0) message("rep ", i, "/", n_reps)
  sim_one()
}))

cat("\n================ RCT validation (", nrow(results), "reps, n =", n_obs, ") ================\n")
cat(sprintf("Modifier detection rate (rule contains W1): %.1f%%\n", 100 * mean(results$detected)))
summ <- function(est, truth, cov) {
  c(bias = mean(est) - truth, sd = sd(est),
    rmse = sqrt(mean((est - truth)^2)), coverage = mean(cov))
}
tab <- rbind(
  `Region V (truth 2)`   = summ(results$v_est,  true_v,   results$v_cov),
  `Region V^c (truth 0)` = summ(results$vc_est, true_vc,  results$vc_cov),
  `Contrast (truth 2)`   = summ(results$con_est, true_con, results$con_cov)
)
print(round(tab, 4))
cat("\n(Coverage near 0.95 and small bias indicate the RCT estimator is valid.)\n")

# ---------------------------------------------------------------------------
# Reference results (100 reps, n = 1000, AIPW region estimator):
#
#   Modifier detection rate (rule contains W1): 100.0%
#                           bias     sd   rmse coverage
#   Region V (truth 2)   -0.0103 0.1025 0.1025     0.92
#   Region V^c (truth 0) -0.0007 0.0845 0.0840     0.98
#   Contrast (truth 2)   -0.0096 0.1352 0.1349     0.94
#
# (An earlier plug-in region estimator gave large bias and ~0.03-0.18 coverage;
#  switching to the doubly-robust AIPW pseudo-outcome restored nominal behaviour.)
# ---------------------------------------------------------------------------
