make_fast_shift_learner <- function() {
  sl3::make_learner(sl3::Stack, c(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new()))
}

make_shift_data <- function(n = 120, n_exposures = 1) {
  W1 <- rbinom(n, 1, 0.5)
  W2 <- rnorm(n)
  A1 <- rnorm(n, mean = 0.3 * W1)
  A2 <- rnorm(n, mean = -0.2 * W2)
  Y <- A1 + 0.5 * W1 + if (n_exposures > 1) 0.2 * A2 else 0 + rnorm(n)

  a <- if (n_exposures == 1) {
    data.frame(A1 = A1)
  } else {
    data.frame(A1 = A1, A2 = A2)
  }

  list(
    w = data.frame(W1 = W1, W2 = W2),
    a = a,
    y = Y
  )
}

test_that("stochastic-shift workflow handles a single exposure with no seed", {
  set.seed(1)
  dat <- make_shift_data(n_exposures = 1)
  learner <- make_fast_shift_learner()

  res <- EffectXshift(
    w = dat$w,
    a = dat$a,
    y = dat$y,
    deltas = 0.1,
    n_folds = 2,
    parallel = FALSE,
    min_obs = 10,
    density_classification = TRUE,
    mu_learner = learner,
    g_learner = learner
  )

  expect_named(res, c(
    "Effect Modification K-Fold Results",
    "Effect Modification Region V Pooled Results",
    "Effect Modification Region V^c Pooled Results",
    "Marginal Shift Results",
    "Selection Diagnostics",
    "Positivity Diagnostics"
  ))
  expect_true(all(c("V", "Vc") %in% res$`Effect Modification K-Fold Results`$RegionType))

  selection <- diagnose_selection(res)
  positivity <- diagnose_positivity(res)

  expect_named(selection, c("fold_level", "exposure_frequency", "rule_frequency", "stability"))
  expect_true(all(c("Hn_Shift_Max", "Flag_Hn_Shift_Max_GT_Threshold") %in% names(positivity)))
  expect_equal(nrow(positivity), 2)
})

test_that("unnamed delta vectors are assigned to exposures in order", {
  set.seed(2)
  dat <- make_shift_data(n_exposures = 2)
  learner <- make_fast_shift_learner()

  res <- EffectXshift(
    w = dat$w,
    a = dat$a,
    y = dat$y,
    deltas = c(0.1, 0.2),
    n_folds = 2,
    parallel = FALSE,
    seed = 10,
    min_obs = 10,
    density_classification = TRUE,
    mu_learner = learner,
    g_learner = learner
  )

  marginal <- res$`Marginal Shift Results`
  expect_setequal(marginal$Variables, c("A1", "A2"))
  expect_setequal(res$`Positivity Diagnostics`$Exposure, c("A1", "A2"))
})

test_that("rct workflow rejects multiple exposure columns", {
  dat <- make_shift_data(n_exposures = 2)

  expect_error(
    EffectXshift(
      w = dat$w,
      a = dat$a,
      y = dat$y,
      deltas = c(0.1, 0.2),
      rct = TRUE,
      parallel = FALSE
    ),
    "exactly one binary exposure"
  )
})

test_that("rct workflow validates known allocation probability", {
  dat <- make_shift_data(n_exposures = 1)

  expect_error(
    EffectXshift(
      w = dat$w,
      a = dat$a,
      y = dat$y,
      deltas = 0,
      rct = TRUE,
      alpha = 1,
      parallel = FALSE
    ),
    "alpha"
  )
})

test_that("tmle_exposhift defaults are callable", {
  set.seed(3)
  data_internal <- data.frame(y = rnorm(10))
  Qn_scaled <- data.frame(noshift = rep(0.45, 10), upshift = rep(0.55, 10))
  Hn <- data.frame(noshift = rep(1, 10), shift = rep(1, 10))

  fit <- tmle_exposhift(
    data_internal = data_internal,
    delta = 0.1,
    Qn_scaled = Qn_scaled,
    Qn_unscaled = Qn_scaled,
    Hn = Hn,
    y = data_internal$y
  )

  expect_true(is.numeric(fit$psi))
})

test_that("calculatePooledEstimate returns appended results, not a function", {
  results <- data.frame(
    Condition = "A1",
    Psi = c(0.5, 0.8),
    Variance = c(0.05, 0.04),
    SE = c(0.22, 0.2),
    `Lower CI` = c(0.1, 0.4),
    `Upper CI` = c(0.9, 1.2),
    `P-value` = c(0.04, 0.01),
    Fold = c("Fold 1", "Pooled TMLE"),
    Type = "Indiv Shift",
    Variables = "A1",
    N = 100,
    Delta = 0.1,
    check.names = FALSE
  )

  pooled <- calculatePooledEstimate(results, n_folds = 3)

  expect_s3_class(pooled, "data.frame")
  expect_false(is.function(pooled))
  expect_true("Inverse Variance Pooled" %in% pooled$Fold)
})
