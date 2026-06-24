# ---------------------------------------------------------------------------
# Breadth validation for the RCT effect-modification workflow:
#   (a) sample-size sweep with a BINARY modifier (true contrast = 2)
#   (b) a CONTINUOUS modifier (effect switches on above a threshold)
# Reports detection, bias, and V - V^c contrast coverage in each setting.
# ---------------------------------------------------------------------------

suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
  library(sl3)
})

run_setting <- function(dgp, n_obs, n_reps, true_contrast, modifier_name, seed0) {
  set.seed(seed0)
  one <- function() {
    d <- dgp(n_obs)
    res <- tryCatch(EffectXshift(
      w = d$w, a = d$a, y = d$y, deltas = 0.1,
      mu_learner = list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new()),
      g_learner = list(Lrnr_ranger$new(num.trees = 100), Lrnr_glm$new()),
      n_folds = 2, parallel = FALSE, seed = 1,
      rct = TRUE, rct_type = "ate", min_obs = 40, max_depth = 1, pval_thresh = 0.2
    ), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    pooled <- res[["Pooled Region Effects"]]; kf <- res[["Effect Modification K-Fold Results"]]
    con <- pooled[pooled$Region == "V - V^c", ]
    data.frame(
      detected = any(grepl(modifier_name, kf$Rule)),
      con_est  = con$Psi,
      cover    = (con$`Lower CI` <= true_contrast & con$`Upper CI` >= true_contrast)
    )
  }
  r <- do.call(rbind, lapply(seq_len(n_reps), function(i) one()))
  data.frame(
    setting = sprintf("%s, n=%d", modifier_name, n_obs),
    reps = nrow(r),
    detection = round(mean(r$detected), 3),
    contrast_mean = round(mean(r$con_est), 3),
    contrast_truth = true_contrast,
    coverage = round(mean(r$cover), 3)
  )
}

# (a) binary modifier W1: effect 2 if W1==1 else 0 -> contrast 2
dgp_binary <- function(n) {
  W1 <- rbinom(n,1,.5); W2 <- rnorm(n); W3 <- rbinom(n,1,.4)
  A <- rbinom(n,1,.5)
  Y <- 1 + 0.5*W2 - 0.3*W3 + ifelse(W1==1,2,0)*A + rnorm(n)
  list(w = data.frame(W1=W1,W2=W2,W3=W3), a = data.frame(A=A), y = Y)
}

# (b) continuous modifier W1: effect 2 if W1>0 else 0 -> among-region contrast ~2
dgp_cont <- function(n) {
  W1 <- rnorm(n); W2 <- rnorm(n); W3 <- rbinom(n,1,.5)
  A <- rbinom(n,1,.5)
  Y <- 1 + 0.5*W2 - 0.3*W3 + ifelse(W1>0,2,0)*A + rnorm(n)
  list(w = data.frame(W1=W1,W2=W2,W3=W3), a = data.frame(A=A), y = Y)
}

reps <- as.integer(Sys.getenv("BREADTH_REPS", "40"))
out <- rbind(
  run_setting(dgp_binary,  500, reps, 2, "W1", 101),
  run_setting(dgp_binary, 1000, reps, 2, "W1", 102),
  run_setting(dgp_binary, 2000, reps, 2, "W1", 103),
  run_setting(dgp_cont,   2000, reps, 2, "W1", 104)
)
cat("\n===== RCT breadth validation =====\n")
print(out, row.names = FALSE)
