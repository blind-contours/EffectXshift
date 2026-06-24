# ===========================================================================
# EffectXshift — RCT subgroup discovery: a self-contained, reproducible demo
# ===========================================================================
# Finds the subgroup of an RCT with the largest treatment-effect modification,
# with honest (sample-split / CV-TMLE) inference on the discovered region.
#
#   install: remotes::install_github("blind-contours/EffectXshift")
# ---------------------------------------------------------------------------

library(EffectXshift)
library(sl3)

set.seed(429)

# ---- Simulate a randomized trial with a real effect-modifying subgroup -----
n  <- 1500
W1 <- rbinom(n, 1, 0.5)        # the true effect modifier (e.g. a biomarker+)
W2 <- rnorm(n)                 # noise covariates
W3 <- rbinom(n, 1, 0.4)
W4 <- rnorm(n)
A  <- rbinom(n, 1, 0.5)        # randomized treatment (alpha = 0.5)

# treatment helps only when W1 == 1 (effect 2); no effect when W1 == 0
tau <- ifelse(W1 == 1, 2, 0)
Y   <- 1 + 0.5 * W2 - 0.3 * W3 + tau * A + rnorm(n)

# ---- Discover the maximally-effect-modified subgroup -----------------------
fit <- EffectXshift(
  w = data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4),
  a = data.frame(A = A),
  y = Y,
  deltas = 0.1,
  rct = TRUE,
  rct_type = "ate",       # subject-level ATE: Q(1,W) - Q(0,W)
  n_folds = 5,
  max_depth = 1,
  seed = 429
)

# Per-fold discovered region rules + held-out effect estimates
print(fit[["Effect Modification K-Fold Results"]])

# Pooled CV-TMLE estimates: effect in V, in V^c, and the V - V^c contrast,
# each with SE, Wald 95% CI, and p-value.
print(fit[["Pooled Region Effects"]])

# Expected: the rule recovers {W1 == 1}; effect ~2 in V, ~0 in V^c, contrast ~2
# with a CI excluding 0 -- the data-adaptively discovered subgroup, validated
# on held-out folds.

# ---- (Optional) single-arm prognostic high-risk region (no contrast) -------
# Yb <- rbinom(n, 1, plogis(-1 + 2 * W1 + 0.5 * W2))
# fit_risk <- EffectXshift(
#   w = data.frame(W1 = W1, W2 = W2, W3 = W3), a = data.frame(A = A), y = Yb,
#   deltas = 0.1, outcome_type = "binary", rct = TRUE, target = "risk",
#   n_folds = 5, max_depth = 1, seed = 429)
# print(fit_risk[["Pooled Region Effects"]])  # held-out region risks (rates)
