#' Summarize Fold-Level Selection Stability
#'
#' @description
#' Summarizes how consistently \code{\link{EffectXshift}} selects exposures and
#' covariate-region rules across folds. Stable selections are easier to interpret
#' than fold-specific rules that frequently change exposure, modifier, threshold,
#' or orientation.
#'
#' @param x An \code{EffectXshift} result object or a data frame like the
#' \code{"Effect Modification K-Fold Results"} table.
#' @param rank Optional integer vector of ranks to summarize. The default
#' \code{NULL} summarizes all ranks.
#'
#' @return A list with fold-level selected regions, exposure frequencies, rule
#' frequencies, and one row per rank summarizing the most common exposure/rule.
#' @export
diagnose_selection <- function(x, rank = NULL) {
  k_fold_results <- extract_kfold_results(x)

  if (!is.data.frame(k_fold_results)) {
    stop("x must be an EffectXshift result object or a k-fold results data frame.", call. = FALSE)
  }
  if (nrow(k_fold_results) == 0) {
    return(list(
      fold_level = data.frame(),
      exposure_frequency = data.frame(),
      rule_frequency = data.frame(),
      stability = data.frame()
    ))
  }

  rule_col <- if ("Modifier" %in% names(k_fold_results)) {
    "Modifier"
  } else if ("Rule" %in% names(k_fold_results)) {
    "Rule"
  } else {
    NA_character_
  }

  region_col <- if ("RegionType" %in% names(k_fold_results)) "RegionType" else NA_character_

  fold_level <- data.frame(
    Fold = k_fold_results$Fold,
    Rank = if ("Rank" %in% names(k_fold_results)) k_fold_results$Rank else 1L,
    Exposure = k_fold_results$Exposure,
    RegionType = if (!is.na(region_col)) k_fold_results[[region_col]] else NA_character_,
    Rule = if (!is.na(rule_col)) k_fold_results[[rule_col]] else NA_character_,
    N_in_Region = if ("N_in_Region" %in% names(k_fold_results)) k_fold_results$N_in_Region else NA_integer_,
    Effect = if ("Effect" %in% names(k_fold_results)) k_fold_results$Effect else NA_real_,
    SE = if ("SE" %in% names(k_fold_results)) k_fold_results$SE else NA_real_,
    stringsAsFactors = FALSE
  )

  if ("RegionType" %in% names(fold_level)) {
    fold_level <- fold_level[is.na(fold_level$RegionType) | fold_level$RegionType == "V", , drop = FALSE]
  }
  if (!is.null(rank)) {
    fold_level <- fold_level[fold_level$Rank %in% rank, , drop = FALSE]
  }

  exposure_frequency <- frequency_table(fold_level, c("Rank", "Exposure"))
  rule_frequency <- frequency_table(fold_level, c("Rank", "Exposure", "Rule"))

  stability <- do.call(
    rbind,
    lapply(split(fold_level, fold_level$Rank), summarize_selection_rank)
  )
  rownames(stability) <- NULL

  list(
    fold_level = fold_level,
    exposure_frequency = exposure_frequency,
    rule_frequency = rule_frequency,
    stability = stability
  )
}

#' Summarize Positivity Diagnostics
#'
#' @description
#' Returns fold/exposure-level clever-covariate diagnostics for the mixed
#' stochastic-shift workflow and adds flags for large estimated density ratios.
#' Large clever covariates indicate weak practical support for the requested
#' shift and should be reviewed before interpreting region estimates.
#'
#' @param x An \code{EffectXshift} result object or a data frame like the
#' \code{"Positivity Diagnostics"} table.
#' @param threshold Optional maximum acceptable value for \code{Hn_Shift_Max}.
#' If \code{NULL}, the stored threshold is used when available, otherwise 10.
#'
#' @return A data frame with one row per fold/exposure and positivity flags.
#' @export
diagnose_positivity <- function(x, threshold = NULL) {
  positivity <- if (is.list(x) && "Positivity Diagnostics" %in% names(x)) {
    x[["Positivity Diagnostics"]]
  } else {
    x
  }

  if (!is.data.frame(positivity)) {
    stop("x must be an EffectXshift result object or a positivity diagnostics data frame.", call. = FALSE)
  }
  if (nrow(positivity) == 0) {
    return(positivity)
  }

  if (is.null(threshold)) {
    threshold <- if ("Threshold" %in% names(positivity)) {
      positivity$Threshold
    } else {
      10
    }
  }
  if (length(threshold) == 1) {
    threshold <- rep(threshold, nrow(positivity))
  }

  positivity$Diagnostic_Threshold <- threshold
  positivity$Flag_Hn_Shift_Max_GT_Threshold <- positivity$Hn_Shift_Max > threshold
  positivity$Flag_Hn_Shift_Q99_GT_Threshold <- positivity$Hn_Shift_Q99 > threshold
  positivity$Flag_Any_Training_Threshold_Violation <-
    positivity$Prop_Hn_Shift_GT_Threshold > 0

  positivity
}

#' Summarize Trial Region Arm Balance
#'
#' @description
#' Summarizes arm counts and observed outcome means in the randomized-trial
#' regions returned by \code{\link{EffectXshift}} when \code{rct = TRUE}. These
#' summaries are descriptive diagnostics for trial interpretation; the adjusted
#' region effects are reported in \code{"Pooled Region Effects"}.
#'
#' @param x An \code{EffectXshift} randomized-trial result object.
#' @param treatment Optional treatment column name. This is only required when
#' \code{x} does not already contain stored \code{"Trial Region Diagnostics"}.
#' @param outcome Outcome column name in the stored region data. Defaults to
#' \code{"y"}, the internal outcome name used by \code{\link{EffectXshift}}.
#'
#' @return A data frame with one row each for \code{V} and \code{V^c}, including
#' total sample size, treatment-arm counts, observed treatment fraction,
#' outcome means by arm, a descriptive unadjusted mean difference, and a small
#' arm-count flag.
#' @export
diagnose_trial_regions <- function(x, treatment = NULL, outcome = "y") {
  if (is.list(x) &&
      "Trial Region Diagnostics" %in% names(x) &&
      is.null(treatment) &&
      identical(outcome, "y")) {
    return(x[["Trial Region Diagnostics"]])
  }

  required <- c("Region V Data", "Region V^c Data")
  if (!is.list(x) || !all(required %in% names(x))) {
    stop(
      "x must be an EffectXshift randomized-trial result object with stored region data.",
      call. = FALSE
    )
  }
  if (is.null(treatment)) {
    stop(
      "treatment must be supplied when Trial Region Diagnostics is not stored in x.",
      call. = FALSE
    )
  }

  make_trial_region_diagnostics(
    region_v_data = x[["Region V Data"]],
    region_vc_data = x[["Region V^c Data"]],
    treatment = treatment,
    outcome = outcome
  )
}

make_trial_region_diagnostics <- function(region_v_data,
                                          region_vc_data,
                                          treatment,
                                          outcome = "y",
                                          small_arm_n = 10) {
  summarize_one_region <- function(data, label) {
    if (!all(c(treatment, outcome) %in% names(data))) {
      stop(
        "Stored region data must contain the treatment and outcome columns.",
        call. = FALSE
      )
    }

    n <- nrow(data)
    if (n == 0) {
      return(data.frame(
        Region = label,
        N = 0L,
        N_Control = 0L,
        N_Treated = 0L,
        N_Missing_Treatment = 0L,
        N_Missing_Outcome = 0L,
        Prop_Treated = NA_real_,
        Outcome_Mean = NA_real_,
        Outcome_Mean_Control = NA_real_,
        Outcome_Mean_Treated = NA_real_,
        Observed_Mean_Diff_Treated_Minus_Control = NA_real_,
        Flag_Small_Arm_N = TRUE,
        stringsAsFactors = FALSE
      ))
    }

    A <- data[[treatment]]
    Y <- data[[outcome]]
    observed_a <- !is.na(A)
    n_control <- sum(A == 0, na.rm = TRUE)
    n_treated <- sum(A == 1, na.rm = TRUE)
    n_observed_a <- sum(observed_a)

    if (n_observed_a > 0 && !all(unique(A[observed_a]) %in% c(0, 1))) {
      warning("treatment should be coded 0/1 for trial region diagnostics.", call. = FALSE)
    }

    mean_or_na <- function(z) {
      if (length(z) == 0 || all(is.na(z))) {
        return(NA_real_)
      }
      mean(z, na.rm = TRUE)
    }

    y_control <- Y[A == 0 & !is.na(A)]
    y_treated <- Y[A == 1 & !is.na(A)]
    mean_control <- mean_or_na(y_control)
    mean_treated <- mean_or_na(y_treated)

    data.frame(
      Region = label,
      N = n,
      N_Control = n_control,
      N_Treated = n_treated,
      N_Missing_Treatment = sum(is.na(A)),
      N_Missing_Outcome = sum(is.na(Y)),
      Prop_Treated = if (n_observed_a > 0) n_treated / n_observed_a else NA_real_,
      Outcome_Mean = mean_or_na(Y),
      Outcome_Mean_Control = mean_control,
      Outcome_Mean_Treated = mean_treated,
      Observed_Mean_Diff_Treated_Minus_Control = mean_treated - mean_control,
      Flag_Small_Arm_N = min(n_control, n_treated) < small_arm_n,
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- rbind(
    summarize_one_region(region_v_data, "V"),
    summarize_one_region(region_vc_data, "V^c")
  )
  rownames(diagnostics) <- NULL
  diagnostics
}

extract_kfold_results <- function(x) {
  if (is.list(x) && "Effect Modification K-Fold Results" %in% names(x)) {
    return(x[["Effect Modification K-Fold Results"]])
  }
  x
}

frequency_table <- function(data, cols) {
  if (nrow(data) == 0) {
    return(data.frame())
  }

  counts <- as.data.frame(table(data[cols], useNA = "ifany"), stringsAsFactors = FALSE)
  names(counts)[names(counts) == "Freq"] <- "N_Folds"
  counts <- counts[counts$N_Folds > 0, , drop = FALSE]
  rownames(counts) <- NULL
  counts[order(-counts$N_Folds), , drop = FALSE]
}

summarize_selection_rank <- function(data) {
  exposure_counts <- sort(table(data$Exposure), decreasing = TRUE)
  rule_counts <- sort(table(data$Rule), decreasing = TRUE)

  data.frame(
    Rank = unique(data$Rank)[1],
    Folds = length(unique(data$Fold)),
    Unique_Exposures = length(unique(data$Exposure)),
    Most_Common_Exposure = names(exposure_counts)[1],
    Most_Common_Exposure_Folds = as.integer(exposure_counts[1]),
    Unique_Rules = length(unique(data$Rule)),
    Most_Common_Rule = names(rule_counts)[1],
    Most_Common_Rule_Folds = as.integer(rule_counts[1]),
    stringsAsFactors = FALSE
  )
}
