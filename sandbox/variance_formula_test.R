# ---------------------------------------------------------------------------
# Pure variance-FORMULA test with ORACLE (true) nuisances.
#
# By plugging in the analytically-known true outcome regression Q, its shifted
# value Q_d, and the true MTP density ratio H, we remove ALL nuisance-estimation
# error. Whatever miscalibration remains is purely the variance FORMULA.
#
# DGP: A1 ~ N(0.3 W2, 1); Y = 1 + A1 + 0.5 A2 + 0.2 A3 + 2 A1 W1 + 0.4 W2 + N(0,1).
#   Shift A1 -> A1 + delta. True per-subject effect Q_d - Q = delta*(1 + 2 W1).
#   Region V = {W1 == 1}: true Psi_V = 3*delta.
#   True density ratio H(a,w) = g(a-delta|w)/g(a|w), g = N(0.3 W2, 1):
#       H = exp(delta*(a - 0.3 W2) - delta^2/2)     [clever covariate at obs a]
#
# We compare three standard errors for the region-V estimate, against the
# empirical SD of the (efficient AIPW) point estimate over many reps:
#   se_package : sd( (H-1)(Y-Q) + (Q_d - mean Q_d) ) / sqrt(n_V)   <- current code
#   se_condEIF : sd( H(Y-Q) + (Q_d - Q) )           / sqrt(n_V)   <- region-conditional EIF
# ---------------------------------------------------------------------------

set.seed(20240615)
B <- 600
n <- 1000
delta <- 1
true_v <- 3 * delta

one <- function() {
  W1 <- rbinom(n, 1, 0.5); W2 <- rnorm(n); W3 <- rnorm(n)
  A1 <- rnorm(n, 0.3 * W2); A2 <- rnorm(n, 0.2 * W3); A3 <- rnorm(n)
  Y  <- 1 + A1 + 0.5 * A2 + 0.2 * A3 + 2 * A1 * W1 + 0.4 * W2 + rnorm(n)

  Q   <- 1 + A1 + 0.5 * A2 + 0.2 * A3 + 2 * A1 * W1 + 0.4 * W2          # true Q(A,W)
  Qd  <- 1 + (A1 + delta) + 0.5 * A2 + 0.2 * A3 + 2 * (A1 + delta) * W1 + 0.4 * W2
  H   <- exp(delta * (A1 - 0.3 * W2) - delta^2 / 2)                     # true density ratio

  idx <- which(W1 == 1)
  Yr <- Y[idx]; Qr <- Q[idx]; Qdr <- Qd[idx]; Hr <- H[idx]; nr <- length(idx)

  # efficient AIPW point estimate of Psi_V (true nuisances)
  psi <- mean((Qdr - Qr) + Hr * (Yr - Qr))

  se_package <- sd((Hr - 1) * (Yr - Qr) + (Qdr - mean(Qdr))) / sqrt(nr)
  se_condEIF <- sd(Hr * (Yr - Qr) + (Qdr - Qr)) / sqrt(nr)

  c(psi = psi, se_package = se_package, se_condEIF = se_condEIF)
}

r <- as.data.frame(t(replicate(B, one())))

emp_sd <- sd(r$psi)
cat(sprintf("\nOracle-nuisance variance-formula test (B=%d, n=%d, region V truth=%.1f)\n", B, n, true_v))
cat(sprintf("  point estimate mean = %.3f   (empirical SD = %.3f)\n", mean(r$psi), emp_sd))
report <- function(name, se) {
  cov <- mean(r$psi - 1.96 * se <= true_v & r$psi + 1.96 * se >= true_v)
  cat(sprintf("  %-12s mean SE = %.3f   SE/empSD = %.2f   coverage = %.2f\n",
              name, mean(se), mean(se) / emp_sd, cov))
}
report("package", r$se_package)
report("condEIF", r$se_condEIF)
cat("\nThe formula whose SE/empSD ~ 1 and coverage ~ 0.95 is the correct one.\n")

# ---------------------------------------------------------------------------
# RESULT (B=600, n=1000):
#   point estimate mean = 2.998   empirical SD = 0.076
#   package   mean SE = 0.160   SE/empSD = 2.09   coverage = 1.00
#   condEIF   mean SE = 0.073   SE/empSD = 0.96   coverage = 0.94
#
# Verdict: with TRUE nuisances the package variance formula is correct (in fact
# ~2x conservative); the region-conditional EIF is exactly efficient. So the
# anti-conservative CIs seen in the full pipeline are NOT a formula bug -- they
# are driven by nuisance-estimation error: the estimator's true SD inflates from
# 0.076 (true nuisances) to ~0.27 (ranger, n=1000), while the IF-based SE stays
# ~0.16. The fix is therefore better/faster nuisances or a resampling-based SE
# that captures nuisance variability -- NOT a change to the variance formula
# (switching to the tighter condEIF would only worsen real-data coverage).
# ---------------------------------------------------------------------------
