# ---------------------------------------------------------------------------
# Validation for the single-arm PROGNOSTIC risk mode (rct = TRUE, target = "risk").
#
# No treatment contrast: the method fits Q(W) = E[Y|W], partitions the smooth
# training risk to discover a HIGH-RISK region, and reports the HELD-OUT
# empirical risk (a prediction, not a causal effect).
#
# DGP: binary outcome with risk driven by W1: P(Y=1|W) = plogis(-1 + 2 W1 + 0.5 W2).
#   High-risk region V = {W1 == 1}; low-risk complement = {W1 == 0}.
# Checks: (1) discovers W1; (2) held-out region-V risk CI covers the true risk;
#         (3) V risk > V^c risk.
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

set.seed(20240623)

# true region risks (large Monte Carlo)
big <- 2e6
W2b <- rnorm(big)
true_v_risk  <- mean(plogis(-1 + 2 * 1 + 0.5 * W2b))  # W1 == 1
true_vc_risk <- mean(plogis(-1 + 2 * 0 + 0.5 * W2b))  # W1 == 0

n_reps <- as.integer(Sys.getenv("RISK_REPS", "40"))
n_obs  <- as.integer(Sys.getenv("RISK_N", "1500"))

sim_one <- function() {
  W1 <- rbinom(n_obs, 1, 0.5); W2 <- rnorm(n_obs); W3 <- rbinom(n_obs, 1, 0.4)
  A  <- rbinom(n_obs, 1, 0.5)  # present but irrelevant for prognostic risk
  Y  <- rbinom(n_obs, 1, plogis(-1 + 2 * W1 + 0.5 * W2))

  res <- tryCatch(EffectXshift(
    w = data.frame(W1 = W1, W2 = W2, W3 = W3), a = data.frame(A = A), y = Y,
    deltas = 0.1, outcome_type = "binary",
    mu_learner = list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new()),
    g_learner = list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new()),
    n_folds = 2, parallel = FALSE, seed = 1,
    rct = TRUE, target = "risk", min_obs = 40, max_depth = 1, pval_thresh = 0.1
  ), error = function(e) conditionMessage(e))
  if (is.character(res)) { message("ERR: ", substr(res,1,60)); return(NULL) }

  pooled <- res[["Pooled Region Effects"]]; kf <- res[["Effect Modification K-Fold Results"]]
  v  <- pooled[pooled$Region == "V", ]; vc <- pooled[pooled$Region == "V^c", ]
  data.frame(
    detected = any(grepl("W1", kf$Rule)),
    v_risk = v$Psi, vc_risk = vc$Psi,
    v_cover = (v$`Lower CI` <= true_v_risk & v$`Upper CI` >= true_v_risk),
    v_gt_vc = (v$Psi > vc$Psi)
  )
}

r <- do.call(rbind, lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0) message("rep ", i, "/", n_reps); sim_one()
}))

cat(sprintf("\n===== RCT risk-mode validation (B=%d, n=%d) =====\n", nrow(r), n_obs))
cat(sprintf("True risks: V (W1=1) = %.3f,  V^c (W1=0) = %.3f\n", true_v_risk, true_vc_risk))
cat(sprintf("Discovered W1 as risk modifier: %.1f%%\n", 100 * mean(r$detected)))
cat(sprintf("V risk est = %.3f (truth %.3f);  V^c risk est = %.3f (truth %.3f)\n",
            mean(r$v_risk), true_v_risk, mean(r$vc_risk), true_vc_risk))
cat(sprintf("Held-out V-risk CI coverage of truth: %.1f%%\n", 100 * mean(r$v_cover)))
cat(sprintf("V risk > V^c risk: %.1f%%\n", 100 * mean(r$v_gt_vc)))
