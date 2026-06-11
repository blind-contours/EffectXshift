library(sl3)

# End-to-end test of the RCT workflow (single binary exposure).
#
# DGP: a randomized binary treatment A (alpha = 0.5) whose effect is modified by
# a binary covariate W1. The treatment effect is ~2 when W1 == 1 and ~0 when
# W1 == 0, so the oracle effect-modification region V should be {W1 == 1} and the
# V - V^c contrast should be clearly positive.
test_that("EffectXshift(rct = TRUE) recovers a binary effect modifier", {
  set.seed(2024)
  n <- 500
  W1 <- rbinom(n, 1, 0.5)
  W2 <- rnorm(n)
  W3 <- rbinom(n, 1, 0.4)
  A <- rbinom(n, 1, 0.5) # randomized treatment

  tau <- ifelse(W1 == 1, 2, 0) # treatment effect modified by W1
  Y <- 1 + 0.5 * W2 - 0.3 * W3 + tau * A + rnorm(n)

  w <- data.frame(W1 = W1, W2 = W2, W3 = W3)
  a <- data.frame(A = A)

  # A learner that can capture the A x W1 interaction (a main-effects GLM cannot).
  mu_learner <- list(Lrnr_ranger$new(num.trees = 200), Lrnr_glm$new())

  res <- EffectXshift(
    w = w, a = a, y = Y,
    deltas = 0.1,
    estimator = "tmle",
    mu_learner = mu_learner,
    g_learner = mu_learner,
    n_folds = 2,
    parallel = FALSE,
    seed = 123,
    rct = TRUE,
    rct_type = "ate",
    min_obs = 30,
    max_depth = 1,
    pval_thresh = 0.2
  )

  # Expected output structure
  expect_named(res, c(
    "Effect Modification K-Fold Results",
    "Pooled Region Effects",
    "Region V Data",
    "Region V^c Data"
  ))

  pooled <- res[["Pooled Region Effects"]]
  expect_true(all(c("V", "V^c", "V - V^c") %in% pooled$Region))
  expect_true(all(c("Psi", "SE", "Lower CI", "Upper CI", "P-value", "N") %in% names(pooled)))

  # The discovered region rule should reference the true modifier W1.
  kfold <- res[["Effect Modification K-Fold Results"]]
  expect_true(any(grepl("W1", kfold$Rule)))

  # Region V (high-effect group) should have a clearly larger effect than V^c,
  # and the oracle contrast should be positive and significant.
  psi_v  <- pooled$Psi[pooled$Region == "V"]
  psi_vc <- pooled$Psi[pooled$Region == "V^c"]
  contrast <- pooled[pooled$Region == "V - V^c", ]

  expect_gt(psi_v, psi_vc)
  expect_gt(contrast$Psi, 0.5)
  expect_lt(contrast$`P-value`, 0.05)
})
