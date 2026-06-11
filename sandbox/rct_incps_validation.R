# ---------------------------------------------------------------------------
# Simulation validation for the incremental-propensity-shift RCT mode
# (EffectXshift(rct = TRUE, rct_type = "incps")).
#
# DGP: same as sandbox/rct_validation.R -- a randomized binary treatment A
# (alpha = 0.5) whose conditional ATE tau(W) = 2 when W1 == 1 and 0 otherwise.
#
# Shifting the treatment probability from alpha to alpha + delta changes the
# mean outcome by delta * E[Q(1,W) - Q(0,W)] = delta * tau(W) per subject.
# Hence the true incremental-shift effect is:
#   Region V  = {W1 == 1}: 2*delta
#   Region V^c = {W1 == 0}: 0
#   Contrast (V - V^c)    : 2*delta
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

set.seed(20240611)

n_reps <- 100
n_obs  <- 1000
delta  <- 0.2
true_v  <- 2 * delta
true_vc <- 0
true_con <- true_v - true_vc

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
      w = w, a = a, y = Y, deltas = delta,
      mu_learner = mu_learner, g_learner = mu_learner,
      n_folds = 2, parallel = FALSE, seed = 1,
      rct = TRUE, rct_type = "incps",
      min_obs = 40, max_depth = 1, pval_thresh = 0.2
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  pooled <- res[["Pooled Region Effects"]]
  kfold  <- res[["Effect Modification K-Fold Results"]]
  get <- function(region) pooled[pooled$Region == region, ]
  v <- get("V"); vc <- get("V^c"); con <- get("V - V^c")

  data.frame(
    detected = any(grepl("W1", kfold$Rule)),
    v_est = v$Psi,  v_cov = (v$`Lower CI`  <= true_v   & v$`Upper CI`  >= true_v),
    vc_est = vc$Psi, vc_cov = (vc$`Lower CI` <= true_vc  & vc$`Upper CI` >= true_vc),
    con_est = con$Psi,
    con_cov = (con$`Lower CI` <= true_con & con$`Upper CI` >= true_con)
  )
}

results <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0) message("rep ", i, "/", n_reps)
  sim_one()
}))

cat("\n===== incps RCT validation (", nrow(results), "reps, n =", n_obs,
    ", delta =", delta, ") =====\n")
cat(sprintf("Modifier detection rate (rule contains W1): %.1f%%\n", 100 * mean(results$detected)))
summ <- function(est, truth, cov) {
  c(bias = mean(est) - truth, sd = sd(est),
    rmse = sqrt(mean((est - truth)^2)), coverage = mean(cov))
}
tab <- rbind(
  `Region V (truth 2*delta)`   = summ(results$v_est,  true_v,   results$v_cov),
  `Region V^c (truth 0)`       = summ(results$vc_est, true_vc,  results$vc_cov),
  `Contrast (truth 2*delta)`   = summ(results$con_est, true_con, results$con_cov)
)
print(round(tab, 4))
cat("\n(Coverage near 0.95 and small bias indicate the incps estimator is valid.)\n")

# ---------------------------------------------------------------------------
# Reference results (100 reps, n = 1000, delta = 0.2, incps = delta * ATE):
#
#   Modifier detection rate (rule contains W1): 100.0%
#                               bias     sd   rmse coverage
#   Region V (truth 2*delta) -0.0021 0.0205 0.0205     0.92
#   Region V^c (truth 0)     -0.0001 0.0169 0.0168     0.98
#   Contrast (truth 2*delta) -0.0019 0.0270 0.0270     0.94
#
# (An earlier implementation differenced the TMLE fluctuation of Q at the
#  observed treatment, which is ~0; it gave 0% detection and ~0 estimates.
#  Recognising that a marginal propensity shift equals delta * ATE fixed it.)
# ---------------------------------------------------------------------------
