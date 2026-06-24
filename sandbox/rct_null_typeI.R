# ---------------------------------------------------------------------------
# NULL / type-I-error check for the RCT effect-modification workflow.
#
# The key credibility question for any subgroup-discovery method: when there is
# NO effect modification, how often does it falsely report a significant
# subgroup? A principled method (discover on training folds, estimate on held-out
# folds via CV-TMLE) should keep this false-discovery rate near the nominal level.
#
# NULL DGP: randomized binary A (alpha = 0.5); a HOMOGENEOUS treatment effect
# (tau = 1 for everyone, no A-by-W interaction) so there is no true effect
# modification. Several covariates (binary + continuous) give the partition many
# chances to find a spurious split; a flexible learner (ranger) can fit spurious
# A:W interactions -> a genuine stress test.
#
# We report the false-discovery rate: the fraction of replicates where the pooled
# V - V^c contrast CI excludes 0 (a falsely "significant" modification subgroup).
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

set.seed(20240623)

n_reps <- as.integer(Sys.getenv("NULL_REPS", "200"))
n_obs  <- as.integer(Sys.getenv("NULL_N", "1000"))

sim_one <- function() {
  W1 <- rbinom(n_obs, 1, 0.5)
  W2 <- rnorm(n_obs)
  W3 <- rbinom(n_obs, 1, 0.4)
  W4 <- rnorm(n_obs)
  W5 <- rbinom(n_obs, 1, 0.5)
  A  <- rbinom(n_obs, 1, 0.5)

  # HOMOGENEOUS effect: tau = 1 everywhere -> NO effect modification.
  # NULL_HARD=TRUE uses a nonlinear baseline so the outcome learner carries
  # realistic estimation noise (more temptation for spurious splits).
  hard <- as.logical(Sys.getenv("NULL_HARD", "FALSE"))
  baseline <- if (hard) {
    0.7 * sin(2 * W2) + 0.5 * W2 * W4 + 0.6 * (W1 == 1) * (W4 > 0) - 0.3 * W3 + 0.4 * W2^2
  } else {
    0.5 * W2 - 0.3 * W3 + 0.4 * W1 + 0.2 * W4
  }
  Y <- 1 + baseline + 1 * A + rnorm(n_obs)

  w <- data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, W5 = W5)
  a <- data.frame(A = A)
  mu_learner <- list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new())

  res <- tryCatch(
    EffectXshift(
      w = w, a = a, y = Y, deltas = 0.1,
      mu_learner = mu_learner, g_learner = mu_learner,
      n_folds = 2, parallel = FALSE, seed = 1,
      rct = TRUE, rct_type = "ate",
      min_obs = 40, max_depth = 1, pval_thresh = 0.1
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(c(claimed = NA, sig = NA))

  pooled <- res[["Pooled Region Effects"]]
  con <- pooled[pooled$Region == "V - V^c", ]
  kf  <- res[["Effect Modification K-Fold Results"]]

  # "claimed a subgroup" = a non-degenerate split was found (a rule referencing a W)
  claimed <- any(nzchar(gsub("NOT\\(|\\)", "", kf$Rule)) & grepl("W", kf$Rule))
  # false positive = significant V - V^c contrast (CI excludes 0)
  sig <- is.finite(con$`Lower CI`) && is.finite(con$`Upper CI`) &&
    (con$`Lower CI` > 0 | con$`Upper CI` < 0)
  c(claimed = as.numeric(claimed), sig = as.numeric(sig))
}

r <- as.data.frame(do.call(rbind, lapply(seq_len(n_reps), function(i) {
  if (i %% 25 == 0) message("rep ", i, "/", n_reps)
  sim_one()
})))
r <- r[!is.na(r$sig), ]

cat(sprintf("\n===== RCT NULL / type-I (B=%d, n=%d, homogeneous effect) =====\n", nrow(r), n_obs))
cat(sprintf("Found a candidate subgroup (any fold):      %.1f%%\n", 100 * mean(r$claimed)))
cat(sprintf("FALSE-POSITIVE rate (V - V^c CI excludes 0): %.1f%%   (target ~5%%)\n",
            100 * mean(r$sig)))
cat("\nA false-positive rate near 5% means the method does not invent subgroups.\n")

# ---------------------------------------------------------------------------
# RESULTS (B=200, n=1000, ranger+glm, 2 folds):
#   Simple baseline (NULL_HARD=FALSE):
#     candidate subgroup found:  0.0%   false-positive rate: 0.0%
#     -> with a learnable homogeneous effect the partition finds no significant
#        difference and never splits.
#   Nonlinear baseline (NULL_HARD=TRUE):
#     candidate subgroup found: 71.0%   false-positive rate: 6.5%  (target ~5%)
#     -> the learner IS tempted to split on noise 71% of the time, but the
#        held-out CV-TMLE inference keeps the false-positive rate near nominal.
#
# Takeaway: the data-adaptive discovery does not inflate false subgroup
# discovery -- sample-splitting + CV-TMLE control the type-I error (~5-6%),
# unlike naive RCT subgroup analysis.
# ---------------------------------------------------------------------------
