#' @title Find Maximum Effect Modifiers (RCT, Single Binary Exposure)
#'
#' @description
#' Estimates effect modification in a randomized trial with a single binary
#' exposure and known (or estimable) randomization probability \code{alpha},
#' under one of two target parameters selected by \code{rct_type}:
#'
#' \itemize{
#'   \item \code{rct_type = "ate"}: the subject-level average treatment effect,
#'         \eqn{Q(1, W_i) - Q(0, W_i)}. The outcome regression is fit with a
#'         Super Learner and updated by a TMLE step using the known propensity.
#'   \item \code{rct_type = "incps"}: an incremental propensity-score shift that
#'         moves the marginal treatment probability from \code{alpha} to
#'         \code{alpha + delta}. For a binary treatment the resulting change in
#'         mean outcome equals \code{delta} times the ATE, so this is estimated
#'         with the same doubly-robust machinery scaled by \code{delta}.
#' }
#'
#' Regions are discovered on the training fold with a recursive partition of the
#' smooth plug-in contrast, and the discovered region \code{V} (plus its true
#' complement \code{V^c}) is estimated on the validation fold using the
#' doubly-robust AIPW pseudo-outcome. Because the randomization probability is
#' known, the AIPW region estimate is unbiased even when the outcome model is
#' misspecified.
#'
#' @param at A training fold \code{data.frame} with columns \code{w_names}, \code{a_name}, \code{outcome}.
#' @param av A validation fold \code{data.frame}.
#' @param delta Numeric scalar; the incremental propensity shift (alpha to alpha + delta) for \code{rct_type = "incps"}.
#' @param a_name Name of the binary exposure (e.g. "A").
#' @param w_names Character vector of baseline covariate names.
#' @param outcome Name of the outcome variable.
#' @param outcome_type One of "continuous", "binary", "count" (for \code{sl3} tasks).
#' @param mu_learner A list of \code{sl3} learners for the outcome regression.
#' @param alpha Known randomization probability; if NULL, estimated as the mean of \code{a_name} in \code{at}.
#' @param top_n Number of top rules to return from the partition search.
#' @param seed Random seed for reproducibility.
#' @param min_obs Minimum number of observations in a valid split branch.
#' @param fold Label for the fold index (for cross-validation).
#' @param max_depth Maximum depth of the partition search tree.
#' @param pval_thresh p-value threshold for accepting a split.
#' @param rct_type Either \code{"ate"} or \code{"incps"}.
#' @param target Either \code{"effect"} (default; the treatment-effect-modification
#'   contrast Q(1,W)-Q(0,W), requires both arms) or \code{"risk"} (a single-arm
#'   PROGNOSTIC risk region: partition the Super-Learner risk Q(W)=E[Y|W] on the
#'   training folds and report the HELD-OUT region risk. No control arm or
#'   propensity is needed; for a region mean the influence-function SE equals the
#'   binomial SE, so Q(W) is used only to discover the region. \code{"risk"}
#'   estimates a prediction, not a causal effect.
#'
#' @return A list with:
#' \item{K_fold_EM_results}{Data frame with 2 rows per discovered region: region V and its complement V^c.}
#' \item{av_q_estimates}{Per-subject AIPW pseudo-outcome on the validation set.}
#' \item{av_hn_estimates}{Per-subject influence-curve contribution on the validation set.}
#' \item{q_region_v, q_region_vc}{Validation pseudo-outcomes in region V and its complement.}
#' \item{g_region_v, g_region_vc}{Validation influence-curve contributions in region V and its complement.}
#' \item{data_region_v, data_region_vc}{The validation rows in V and V^c.}
#' \item{data}{The full validation set with appended estimation columns.}
#'
#' @examples
#' \dontrun{
#' n <- 500
#' W1 <- rbinom(n, 1, 0.5); W2 <- rnorm(n)
#' A  <- rbinom(n, 1, 0.5)
#' Y  <- 1 + ifelse(W1 == 1, 2, 0) * A + 0.5 * W2 + rnorm(n)
#' dat <- data.frame(W1 = W1, W2 = W2, A = A, y = Y)
#' res <- find_max_effect_mods_rct(
#'   at = dat, av = dat, delta = 0.1, a_name = "A",
#'   w_names = c("W1", "W2"), outcome = "y", outcome_type = "continuous",
#'   mu_learner = list(sl3::Lrnr_ranger$new(), sl3::Lrnr_glm$new()),
#'   seed = 1, min_obs = 30, fold = 1, max_depth = 1, rct_type = "ate"
#' )
#' }
#'
#' @export
find_max_effect_mods_rct <- function(
    at,
    av,
    delta,
    a_name,
    w_names,
    outcome,
    outcome_type,
    mu_learner,
    alpha = NULL,
    top_n = 3,
    seed,
    min_obs,
    fold,
    max_depth = 2,
    pval_thresh = 0.05,
    rct_type = c("ate","incps"),
    target = c("effect", "risk")
) {
  rct_type <- match.arg(rct_type)
  target   <- match.arg(target)
  future::plan(future::sequential, gc = TRUE)
  set.seed(seed)

  # ---------------------------------------------------------
  # 0) Simple mean-imputation for demonstration
  # ---------------------------------------------------------
  at <- at %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(w_names),
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    ))
  av <- av %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(w_names),
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    ))

  # ============================================================
  # PER-SUBJECT TARGET QUANTITY
  #   target = "effect": contrast Q(1,W) - Q(0,W)   (needs both arms)
  #   target = "risk"  : single-arm risk Q(W)=E[Y|W] (prognostic; no control arm)
  # ============================================================
  if (target == "effect") {

  # ---------------------------------------------------------
  # 1) alpha if not given
  # ---------------------------------------------------------
  if (is.null(alpha)) {
    alpha <- mean(at[[a_name]])
  }
  # Guardrail: a contrast needs a usable control arm. Warn when the smaller arm
  # has too few EVENTS to anchor Q(0,W) -- the discovered regions will be driven
  # by between-fold noise rather than real effect modification.
  .min_arm_events <- min(
    sum(at[[outcome]][at[[a_name]] == 1], na.rm = TRUE),
    sum(at[[outcome]][at[[a_name]] == 0], na.rm = TRUE)
  )
  if (is.finite(.min_arm_events) && .min_arm_events < 10) {
    warning(sprintf(
      paste0("find_max_effect_mods_rct: only %.0f events in the smaller treatment ",
             "arm (fold %s). The contrast Q(1,W)-Q(0,W) is thinly anchored and ",
             "regions may reflect fold-to-fold noise. Consider target = 'risk' for a ",
             "single-arm prognostic analysis."),
      .min_arm_events, as.character(fold)))
  }
  # delta only shifts the propensity in the incremental-propensity-shift mode;
  # the ATE mode ignores it, so only enforce the [0, 1] bound for "incps".
  if (rct_type == "incps" && ((alpha + delta) < 0 || (alpha + delta) > 1)) {
    stop("alpha + delta must be in [0,1], got: ", alpha + delta)
  }

  # ---------------------------------------------------------
  # 2) Fit outcome regression
  # ---------------------------------------------------------
  covars <- c(a_name, w_names)
  train_task <- sl3::sl3_Task$new(
    data = at,
    covariates = covars,
    outcome = outcome,
    outcome_type = outcome_type
  )
  sl_mu <- sl3::Lrnr_sl$new(
    learners = mu_learner,
    metalearner = sl3::Lrnr_nnls$new()
  )
  mu_fit <- sl_mu$train(train_task)

  # Generate separate predictions Q(1,W) and Q(0,W) for training
  at_for_1 <- at
  at_for_1[[a_name]] <- 1
  at_for_0 <- at
  at_for_0[[a_name]] <- 0

  train_task_A1 <- sl3::sl3_Task$new(
    data = at_for_1,
    covariates = covars,
    outcome = outcome,
    outcome_type = outcome_type
  )
  train_task_A0 <- sl3::sl3_Task$new(
    data = at_for_0,
    covariates = covars,
    outcome = outcome,
    outcome_type = outcome_type
  )

  Q1_init_train <- mu_fit$predict(train_task_A1)
  Q0_init_train <- mu_fit$predict(train_task_A0)

  at$Q1_init <- Q1_init_train
  at$Q0_init <- Q0_init_train

  # Now predictions Q(1,W), Q(0,W) for validation
  av_for_1 <- av
  av_for_1[[a_name]] <- 1
  av_for_0 <- av
  av_for_0[[a_name]] <- 0

  valid_task_A1 <- sl3::sl3_Task$new(
    data = av_for_1,
    covariates = covars,
    outcome = outcome,
    outcome_type = outcome_type
  )
  valid_task_A0 <- sl3::sl3_Task$new(
    data = av_for_0,
    covariates = covars,
    outcome = outcome,
    outcome_type = outcome_type
  )

  Q1_init_valid <- mu_fit$predict(valid_task_A1)
  Q0_init_valid <- mu_fit$predict(valid_task_A0)

  av$Q1_init <- Q1_init_valid
  av$Q0_init <- Q0_init_valid

  # ---------------------------------------------------------
  # Helper functions (TMLE)
  # ---------------------------------------------------------
  get_fluctuation_model <- function(Y, offset_term, clever_cov, data, outcome_type) {
    if (outcome_type %in% c("binary", "binomial")) {
      # Bound Q away from {0,1}: flexible learners can predict exactly 0/1 for a
      # rare outcome, and qlogis(0|1) = +/-Inf would crash the fluctuation glm.
      offset_term <- pmin(pmax(offset_term, 1e-6), 1 - 1e-6)
      glmfit <- stats::glm(
        formula = as.formula(
          paste(Y, "~ -1 + offset(qlogis(offset_term)) +", clever_cov)
        ),
        family = binomial(),
        data = data
      )
      eps <- coef(glmfit)[clever_cov]
      if (is.na(eps)) eps <- 0
      Q_star <- plogis(qlogis(offset_term) + eps * data[[clever_cov]])
    } else {
      glmfit <- stats::glm(
        formula = as.formula(
          paste(Y, "~ -1 + offset(offset_term) +", clever_cov)
        ),
        family = gaussian(),
        data = data
      )
      eps <- coef(glmfit)[clever_cov]
      if (is.na(eps)) eps <- 0
      Q_star <- offset_term + eps * data[[clever_cov]]
    }
    list(Q_star = Q_star, eps = eps)
  }

  tmle_ate_binary <- function(data, outcome, a_name, alpha, outcome_type) {
    A <- data[[a_name]]
    Y <- data[[outcome]]

    # offset depends on observed A
    data$Q_init <- ifelse(A == 1, data$Q1_init, data$Q0_init)

    # Clever covariate for ATE
    data$H_ATE <- (A / alpha) - ((1 - A)/(1 - alpha))

    fluc <- get_fluctuation_model(
      Y = outcome,
      offset_term = data$Q_init,
      clever_cov = "H_ATE",
      data = data,
      outcome_type = outcome_type
    )
    eps <- fluc$eps
    Q_star <- fluc$Q_star

    # Now get CF predictions
    data_1 <- data
    data_1[[a_name]] <- 1
    data_1$offset_cf <- data_1$Q1_init
    data_1$H_ATE_cf <- 1 / alpha

    data_0 <- data
    data_0[[a_name]] <- 0
    data_0$offset_cf <- data_0$Q0_init
    data_0$H_ATE_cf <- -1 / (1 - alpha)

    if (outcome_type %in% c("binary", "binomial")) {
      bcf1 <- pmin(pmax(data_1$offset_cf, 1e-6), 1 - 1e-6)
      bcf0 <- pmin(pmax(data_0$offset_cf, 1e-6), 1 - 1e-6)
      data_1$Q_star_1 <- plogis(qlogis(bcf1) + eps * data_1$H_ATE_cf)
      data_0$Q_star_0 <- plogis(qlogis(bcf0) + eps * data_0$H_ATE_cf)
    } else {
      data_1$Q_star_1 <- data_1$offset_cf + eps * data_1$H_ATE_cf
      data_0$Q_star_0 <- data_0$offset_cf + eps * data_0$H_ATE_cf
    }

    effect_i <- data_1$Q_star_1 - data_0$Q_star_0
    psi_hat <- mean(effect_i, na.rm = TRUE)
    partial_eif <- data$H_ATE * (Y - Q_star)
    full_eif <- partial_eif + (effect_i - psi_hat)

    list(
      effect = effect_i,
      eif = full_eif,
      psi_hat = psi_hat
    )
  }

  # ---------------------------------------------------------
  # 3) Depending on rct_type, define "Effect" and "IC"
  # ---------------------------------------------------------
  # For each subject we keep two quantities:
  #   * Effect: the smooth plug-in contrast (Q(1,W) - Q(0,W), or the shift
  #     analogue). Used to *discover* effect-modification regions on the
  #     training fold -- it is stable and well suited to recursive partitioning.
  #   * Dstar: the doubly-robust AIPW pseudo-outcome
  #       Effect + (clever covariate) * (Y - Q),
  #     i.e. the influence function recentred at the point estimate. Used to
  #     *estimate* region effects on the validation fold. Because the
  #     randomization probability is known, Dstar is unbiased for the region
  #     effect even when the outcome model is misspecified; the plug-in contrast
  #     alone inherits the bias of the outcome learner.
  if (rct_type == "ate") {
    # A) ATE
    at_res <- tmle_ate_binary(at, outcome, a_name, alpha, outcome_type)
    av_res <- tmle_ate_binary(av, outcome, a_name, alpha, outcome_type)

    at$Effect <- at_res$effect
    at$IC_full <- at_res$eif
    at$Dstar  <- at_res$eif + at_res$psi_hat

    av$Effect <- av_res$effect
    av$IC_full <- av_res$eif
    av$Dstar  <- av_res$eif + av_res$psi_hat
  } else {
    # B) Incremental propensity shift (alpha -> alpha + delta).
    # When a binary treatment's marginal assignment probability is shifted from
    # alpha to alpha + delta, the resulting change in mean outcome is
    #   (alpha + delta) E[Q(1,W)] + (1 - alpha - delta) E[Q(0,W)]
    #     - { alpha E[Q(1,W)] + (1 - alpha) E[Q(0,W)] }
    #   = delta * E[Q(1,W) - Q(0,W)],
    # i.e. delta times the ATE. Per subject the effect is delta * (Q(1,W_i) -
    # Q(0,W_i)). We therefore reuse the (doubly-robust) ATE machinery and scale
    # every quantity by delta.
    at_res <- tmle_ate_binary(at, outcome, a_name, alpha, outcome_type)
    av_res <- tmle_ate_binary(av, outcome, a_name, alpha, outcome_type)

    at$Effect  <- delta * at_res$effect
    at$IC_full <- delta * at_res$eif
    at$Dstar   <- delta * (at_res$eif + at_res$psi_hat)

    av$Effect  <- delta * av_res$effect
    av$IC_full <- delta * av_res$eif
    av$Dstar   <- delta * (av_res$eif + av_res$psi_hat)
  }

  } else {
    # =========================================================
    # target = "risk": single-arm prognostic risk Q(W) = E[Y|W]
    # ---------------------------------------------------------
    # No treatment node is intervened on, so there is no propensity term and no
    # control arm is required. We fit Q on the covariates, partition the smooth
    # TRAINING risk score to discover the high-risk region, and report the
    # HELD-OUT empirical risk on validation. For an (unconditional) region mean
    # the efficient influence function is 1(W in V)/pi_V * (Y - psi_V), so the
    # IF-based SE equals the binomial SE -- Q(W) buys DISCOVERY, not efficiency.
    # =========================================================
    q_task_tr <- sl3::sl3_Task$new(
      data = at, covariates = w_names, outcome = outcome, outcome_type = outcome_type)
    sl_q <- sl3::Lrnr_sl$new(learners = mu_learner, metalearner = sl3::Lrnr_nnls$new())
    q_fit <- sl_q$train(q_task_tr)
    at$Qhat <- q_fit$predict(q_task_tr)
    av$Qhat <- q_fit$predict(sl3::sl3_Task$new(
      data = av, covariates = w_names, outcome = outcome, outcome_type = outcome_type))

    # training: partition the smooth risk score (Effect = Q-hat);
    # validation: held-out empirical outcome (Effect = Y) -> honest region rate.
    at$Effect <- at$Qhat;            at$IC_full <- at[[outcome]]; at$Dstar <- at[[outcome]]
    av$Effect <- av[[outcome]];      av$IC_full <- av[[outcome]]; av$Dstar <- av[[outcome]]
  }

  # ---------------------------------------------------------
  # 4) Partition (Training set)
  # ---------------------------------------------------------
  at_df <- data.frame(at)
  at_df$y <- at_df$Effect
  at_df$ic <- at_df$IC_full

  calculate_p_value <- function(left_eff, left_ic, right_eff, right_ic) {
    left_mean <- mean(left_eff, na.rm = TRUE)
    right_mean <- mean(right_eff, na.rm = TRUE)

    left_se <- sqrt(var(left_ic, na.rm = TRUE)/length(left_ic))
    right_se <- sqrt(var(right_ic, na.rm = TRUE)/length(right_ic))

    t_stat <- (left_mean - right_mean) / sqrt(left_se^2 + right_se^2)

    df <- min(length(left_ic), length(right_ic)) - 1
    if (df > 30) {
      p_value <- 2 * pnorm(-abs(t_stat))
    } else {
      p_value <- 2 * pt(-abs(t_stat), df = df)
    }
    p_value
  }

  split_data <- function(data, split_var, split_point) {
    unique_vals <- sort(unique(data[[split_var]]))
    is_binary <- (length(unique_vals) == 2 && all(unique_vals == c(0, 1)))

    if (is_binary) {
      left <- subset(data, data[[split_var]] == 0)
      right <- subset(data, data[[split_var]] == 1)
    } else {
      left <- subset(data, data[[split_var]] <= split_point)
      right <- subset(data, data[[split_var]] > split_point)
    }
    list(left = left, right = right)
  }

  find_best_split <- function(
    data, split_vars, outcome = "y",
    min_obs = 10, pval_thresh = 0.05
  ) {
    best_split <- NULL
    best_abs_diff <- 0

    for (var in split_vars) {
      vals <- sort(unique(data[[var]]))
      is_binary <- (length(vals) == 2 && all(vals == c(0, 1)))

      for (sp in vals) {
        if (is_binary) {
          data$split_indicator <- ifelse(data[[var]] == sp, 1, 0)
        } else {
          data$split_indicator <- ifelse(data[[var]] <= sp, 1, 0)
        }

        n_left <- sum(data$split_indicator == 1)
        n_right <- sum(data$split_indicator == 0)
        if (n_left < min_obs || n_right < min_obs) next

        left_eff <- data$y[data$split_indicator == 1]
        right_eff <- data$y[data$split_indicator == 0]
        left_ic <- data$ic[data$split_indicator == 1]
        right_ic <- data$ic[data$split_indicator == 0]

        pval <- calculate_p_value(left_eff, left_ic, right_eff, right_ic)
        diff <- mean(left_eff, na.rm = TRUE) - mean(right_eff, na.rm = TRUE)

        if (pval <= pval_thresh) {
          if (abs(diff) > best_abs_diff) {
            best_abs_diff <- abs(diff)
            best_split <- list(
              variable = var, point = sp, p_value = pval, diff = diff
            )
          }
        }
      }
    }
    data$split_indicator <- NULL
    best_split
  }

  recursive_split_all_rules <- function(
    data, split_vars, depth = 0, max_depth = 1, outcome = "y",
    path = "", min_obs = 10, pval_thresh = 0.05
  ) {
    if (depth == max_depth || nrow(data) == 0) {
      region_mean <- mean(data[[outcome]], na.rm = TRUE)
      return(list(list(
        RegionMean = region_mean,
        N = nrow(data),
        Rule = path,
        Depth = depth
      )))
    }

    best_spl <- find_best_split(data, split_vars, outcome, min_obs, pval_thresh)
    if (is.null(best_spl)) {
      region_mean <- mean(data[[outcome]], na.rm = TRUE)
      return(list(list(
        RegionMean = region_mean,
        N = nrow(data),
        Rule = path,
        Depth = depth
      )))
    }

    sp_data <- split_data(data, best_spl$variable, best_spl$point)
    is_binary <- (length(unique(data[[best_spl$variable]])) == 2)

    if (is_binary) {
      left_rule <- paste0(path, best_spl$variable, "==", best_spl$point, " & ")
      right_rule <- paste0(path, best_spl$variable, "!=", best_spl$point, " & ")
    } else {
      left_rule <- paste0(path, best_spl$variable, "<=", best_spl$point, " & ")
      right_rule <- paste0(path, best_spl$variable, ">", best_spl$point, " & ")
    }

    left_res <- recursive_split_all_rules(
      sp_data$left, split_vars, depth + 1, max_depth, outcome,
      left_rule, min_obs, pval_thresh
    )
    right_res <- recursive_split_all_rules(
      sp_data$right, split_vars, depth + 1, max_depth, outcome,
      right_rule, min_obs, pval_thresh
    )

    c(left_res, right_res)
  }

  # run the partition
  tree_results <- recursive_split_all_rules(
    data = at_df,
    split_vars = w_names,
    max_depth = max_depth,
    outcome = "y",
    min_obs = min_obs,
    pval_thresh = pval_thresh
  )
  iteration_df <- as.data.frame(do.call(rbind, tree_results))
  iteration_df$Exposure <- a_name
  iteration_df$RegionMean <- as.numeric(iteration_df$RegionMean)
  iteration_df$N <- as.integer(iteration_df$N)

  # Region V is the highest-effect leaf, with V^c its complement. Ordering by
  # the (signed) effect rather than absolute deviation gives a consistent region
  # orientation across folds, so pooling the validation folds does not mix
  # high- and low-effect subgroups.
  if (nrow(iteration_df) > 1) {
    idx <- order(iteration_df$RegionMean, decreasing = TRUE)
    selected_items <- iteration_df[idx[1:min(top_n, length(idx))], ]
  } else {
    selected_items <- iteration_df
  }

  # Now evaluate discovered region(s) in validation. Region effects and their
  # influence-curve-based variances are computed from the doubly-robust AIPW
  # pseudo-outcome (Dstar), not the plug-in contrast, so the point estimate and
  # its standard error are consistent with one another.
  av_df <- data.frame(av)
  av_df$Effect <- av$Dstar
  av_df$ic <- av$Dstar

  # We'll build up a data frame with 2 rows for each discovered region: region vs complement
  results_rows <- list()
  row_counter <- 1

  for (i in seq_len(nrow(selected_items))) {
    item <- selected_items[i, ]
    rule <- item$Rule
    if (nchar(rule) > 0) {
      r_clean <- sub("&\\s*$", "", rule)
      av_df$Indicator <- with(av_df, eval(parse(text = paste0("(", r_clean, ")"))))
      av_df$Indicator <- ifelse(is.na(av_df$Indicator), FALSE, av_df$Indicator)
    } else {
      # if there's no rule, entire set is region
      av_df$Indicator <- TRUE
    }

    sub_idx <- which(av_df$Indicator == TRUE)
    eff_sub <- av_df$Effect[sub_idx]
    ic_sub <- av_df$ic[sub_idx]

    effect_est <- mean(eff_sub, na.rm = TRUE)
    var_est <- var(ic_sub, na.rm = TRUE) / length(ic_sub)
    se_est <- sqrt(var_est)
    ci_lower <- effect_est - 1.96 * se_est
    ci_upper <- effect_est + 1.96 * se_est

    # row for the discovered region V
    region_df <- data.frame(
      Exposure = a_name,
      RegionType = "V",                   # region
      Rule = rule,
      Effect = effect_est,
      SE = se_est,
      LowerCI = ci_lower,
      UpperCI = ci_upper,
      Fold = fold,
      Rank = i,
      N_in_Region = length(sub_idx)
    )

    colnames(region_df)[3] <- "Rule"
    results_rows[[row_counter]] <- region_df
    row_counter <- row_counter + 1

    # row for the complement V^c
    complement_idx <- which(av_df$Indicator == FALSE)
    eff_sub_c <- av_df$Effect[complement_idx]
    ic_sub_c <- av_df$ic[complement_idx]

    if (length(complement_idx) > 0) {
      effect_est_c <- mean(eff_sub_c, na.rm = TRUE)
      var_est_c <- var(ic_sub_c, na.rm = TRUE) / length(ic_sub_c)
      se_est_c <- sqrt(var_est_c)
      ci_lower_c <- effect_est_c - 1.96 * se_est_c
      ci_upper_c <- effect_est_c + 1.96 * se_est_c

      complement_df <- data.frame(
        Exposure = a_name,
        RegionType = "Vc",                 # complement
        Rule = paste0("NOT(", rule, ")"),  # indicate complement
        Effect = effect_est_c,
        SE = se_est_c,
        LowerCI = ci_lower_c,
        UpperCI = ci_upper_c,
        Fold = fold,
        Rank = i,
        N_in_Region = length(complement_idx)
      )
      results_rows[[row_counter]] <- complement_df
      row_counter <- row_counter + 1
    }
  }

  # rbind them
  K_fold_EM_results <- do.call(rbind, results_rows)

  # for top_n=1, define "main" region from selected_items[1, ]
  # build the subset for V and V^c for the "first" discovered region
  if (nrow(selected_items) > 0) {
    main_rule <- selected_items$Rule[1]
    main_rule_clean <- sub("&\\s*$", "", main_rule)

    if (nchar(main_rule_clean) > 0) {
      av_df$Indicator <- with(av_df, eval(parse(text = paste0("(", main_rule_clean, ")"))))
      av_df$Indicator <- ifelse(is.na(av_df$Indicator), FALSE, av_df$Indicator)
    } else {
      av_df$Indicator <- TRUE
    }
  } else {
    av_df$Indicator <- TRUE
  }

  region_data_v <- av_df[av_df$Indicator == TRUE, , drop = FALSE]
  region_data_vc <- av_df[av_df$Indicator == FALSE, , drop = FALSE]

  g_region_v <- region_data_v$ic
  g_region_vc <- region_data_vc$ic
  q_region_v <- region_data_v$Effect
  q_region_vc <- region_data_vc$Effect

  # final
  return(list(
    K_fold_EM_results = K_fold_EM_results,
    av_q_estimates = av_df$Effect,
    av_hn_estimates = av_df$ic,
    g_region_v = g_region_v,
    g_region_vc = g_region_vc,
    q_region_v = q_region_v,
    q_region_vc = q_region_vc,
    data_region_v = region_data_v,
    data_region_vc = region_data_vc,
    data = av_df
  ))
}
