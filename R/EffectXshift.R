#' @title EffectXshift
#'
#' @description This packages aims to find, in a mixed exposure with many covariates and a future outcome, the
#' exposure-covariate region combination that maximizes the differential impact of stochastic shift intervention. In a
#' training fold we find the region in the covariate space that maximizes the average difference in stochastic shift interventions,
#' compared between the regions. This outputs an exposure-covariate region pairing. Then in an estimation sample, we estimate
#' a stochastic shift intervention using targeted learning in each level of the covariate, shifting the discovered exposure.
#' This is done in a CV-TMLE procedure where each fold is used as validation and the complementary folds are used as training.
#' This package outputs the targeted estimates for the discovered exposure in the levels of the discovered covariate both at the k-fold
#' specific level and pooled across folds, which estimates the overall oracle parameter.
#'
#'
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#' baseline covariates. These variables are measured before exposures.
#' @param a \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple exposures.
#' @param y \code{numeric} vector of observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in exposures to
#' define the target parameter, with respect to the scale of the exposures (A). If adaptive_delta
#' is true, these values will be reduced.
#' @param estimator The type of estimator to fit: \code{"tmle"} for targeted
#' maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param fluctuation Method used in the targeting step for TML estimation: "standard" or "weighted".
#' This determines where to place the auxiliary covariate in the logistic tilting regression.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param pi_learner Learners for fitting Super Learner ensembles to the g-mechanism density estimation
#' g(A|W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param n_folds Number of folds to use in cross-validation, default is 2.
#' @param outcome_type Data type of the outcome, default is "continuous".
#' @param parallel Whether to parallelize across cores (default: TRUE).
#' @param parallel_type Type of parallelization to use if parallel is TRUE:
#' "multi_session" (default), "multicore", or "sequential".
#' @param num_cores Number of CPU cores to use in parallelization (default: 2).
#' @param seed \code{numeric} seed value to be passed to all functions.
#' @param hn_trunc_thresh Truncation level for the clever covariate (default: 10).
#' @param adaptive_delta If TRUE, reduces the user-specified delta until
#' @param min_obs Minimum number of observations allowed in a covariate region
#' the Hn calculated for a shift does not have any observation greater
#' than hn_trunc_thresh (default: FALSE).
#' @param top_n top number of effect modifier-exposure pairs to estimate
#'
#' @return An S3 object of class \code{EffectXshift} containing the results of the
#' procedure to compute a TML or one-step estimate of the counterfactual mean
#' under a modified treatment policy that shifts a continuous-valued exposure
#' by a scalar amount \code{delta} in subregions of the exposure space.
#' These exposures are data-adaptively identified using the CV-TMLE procedure.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom stringr str_count
#' @import furrr
#' @importFrom purrr map
#' @importFrom data.table rbindlist

EffectXshift <- function(w,
                         a,
                         y,
                         deltas,
                         estimator = "tmle",
                         fluctuation = "standard",
                         g_learner = NULL,
                         mu_learner = NULL,
                         n_folds = 2,
                         outcome_type = "continuous",
                         parallel = TRUE,
                         parallel_type = "multi_session",
                         num_cores = 2,
                         seed = seed,
                         hn_trunc_thresh = 10,
                         adaptive_delta = FALSE,
                         top_n = 1,
                         min_obs = 20,
                         density_classification = FALSE) {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.data.frame(w)) w <- as.data.frame(w)
  w_names <- colnames(w)
  if (is.null(w_names)) {
    w_names <- paste0("w", seq_len(ncol(w)))
    colnames(w) <- w_names
  }


  # coerce W to matrix and, if no names in W, assign them generically
  a <- data.frame(a)
  a_names <- colnames(a)

  if (is.null(a_names)) {
    a_names <- paste0("a", seq_len(ncol(a)))
    colnames(a) <- a_names
  }

  # If NULL create default learners

  if (is.null(g_learner)) {
    sls <- create_sls()
    g_learner <- sls$g_learner
  }

  if (is.null(mu_learner)) {
    sls <- create_sls()
    mu_learner <- sls$mu_learner
  }

  # Set up parallel type

  if (parallel == TRUE) {
    if (parallel_type == "multi_session") {
      future::plan(future::multisession,
        workers = num_cores,
        gc = TRUE
      )
    } else {
      future::plan(future::multicore,
        workers = num_cores,
        gc = TRUE
      )
    }
  } else {
    future::plan(future::sequential,
      gc = TRUE
    )
  }

  # Create internal data

  data_internal <- data.table::data.table(w, a, y)
  `%notin%` <- Negate(`%in%`)


  # Create folds for CV procedure
  data_internal$folds <- create_cv_folds(n_folds, data_internal$y)

  # Use g-computation and SL to get the individual exposure effects of shifting
  # Then regress these effect vectors onto the covariate space using a greedy partitioning
  # algorithm which finds the region which maximizes the average effect difference in the region vs.
  # out of the region

  # k_fold_results <- list()
  # exposures_shift_g_results <- list()
  # exposures_shift_q_results <- list()
  # g_ests_v_vc <- list()
  # q_ests_v_vc <- list()

  fold_basis_results <- furrr::future_map(unique(data_internal$folds),
    function(fold_k) {
      at <- data_internal[data_internal$folds != fold_k, ]
      av <- data_internal[data_internal$folds == fold_k, ]

      effect_mod_results <- find_max_effect_mods(
        at = at,
        av = av,
        deltas = deltas,
        a_names = a_names,
        w_names = w_names,
        outcome = "y",
        outcome_type = outcome_type,
        mu_learner = mu_learner,
        g_learner = g_learner,
        seed = seed,
        top_n = top_n,
        min_obs = min_obs,
        fold = fold_k,
        density_classification
      )

      k_fold_results <- effect_mod_results$`K-fold_EM_results`
      exposures_shift_q_results <- effect_mod_results$av_q_estimates
      exposures_shift_g_results <- effect_mod_results$av_hn_estimates
      g_region_v <- effect_mod_results$g_region_v
      g_region_vc <- effect_mod_results$g_region_vc
      q_region_v <- effect_mod_results$q_region_v
      q_region_vc <- effect_mod_results$q_region_vc
      data_region_v <- effect_mod_results$region_data_v
      data_region_vc <- effect_mod_results$region_data_vc
      data <- effect_mod_results$data


      results_list <- list(
        k_fold_results,
        exposures_shift_q_results,
        exposures_shift_g_results,
        g_region_v,
        g_region_vc,
        q_region_v,
        q_region_vc,
        data_region_v,
        data_region_vc,
        data
      )

      names(results_list) <- c(
        "k_fold_results",
        "exposure_shift_q",
        "exposure_shift_g",
        "g_estimates_region_v",
        "g_estimates_region_vc",
        "q_estimates_region_v",
        "q_estimates_region_vc",
        "data_region_v",
        "data_region_vc",
        "data"
      )

      results_list

    },
    .options = furrr::furrr_options(seed = seed, packages = "EffectXshift")
  )

  # effect_mod_fold_results <- list()

  # Now go through the results tables for each fold

  # fold_effectXshift_results <- furrr::future_map(
  #   unique(data_internal$folds), function(fold_k) {
  #     ## get the effect modifier and exposure pair with levels
  #
  #     fold_intxn_results <- fold_basis_results[[fold_k]]
  #
  #     ## this gives a table indexed by the exposure that was found with two rows
  #     ## that correspond to the levels of the modifier
  #
  #     for (rank in 1:length(fold_intxn_results)) {
  #
  #       ## get the rank which should be a df with two rows, one for each level
  #       rank_results <- fold_intxn_results[[rank]]
  #
  #       ## for each level of the modifier
  #       for (i in 1:nrow(rank_results)) {
  #         rank_row <- rank_results[i, ]
  #
  #         ## above this is a single row that lists the exposure, rank, and rule
  #         ## for modifier
  #
  #         exposure <- rank_row$Exposure
  #         effect_mod_rule <- unlist(rank_row$Rule)
  #         rank <- rank_row$Rank
  #
  #         ## get name of modifier out of rule for future organization
  #         effect_modifier <- Filter(function(item) grepl(item, effect_mod_rule), w_names)
  #
  #         ## get our training and validation folds
  #         at <- data_internal[data_internal$folds != fold_k, ]
  #         av <- data_internal[data_internal$folds == fold_k, ]
  #
  #
  #         ## get delta for the exposure - how much to shift
  #         delta <- deltas[[exposure]]
  #
  #
  #         ## get lower and upper bound so we don't shift too far
  #         lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
  #         upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))
  #
  #
  #         ## subset to observations that meet the current effect modifier rule
  #         subset_at <- subset(at, eval(parse(text = effect_mod_rule)))
  #         subset_av <- subset(av, eval(parse(text = effect_mod_rule)))
  #
  #
  #         ## we now estimate the density of exposure given covariates in the
  #         ## subpopulation
  #         ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
  #           exposure = exposure,
  #           delta = delta,
  #           g_learner = g_learner,
  #           covars = w_names,
  #           av = subset_av,
  #           at = subset_at,
  #           adaptive_delta = adaptive_delta,
  #           hn_trunc_thresh = hn_trunc_thresh,
  #           use_multinomial = FALSE,
  #           lower_bound = lower_bound,
  #           upper_bound = upper_bound,
  #           outcome_type = "continuous",
  #           density_type = "sl",
  #           n_bins = n_bins,
  #           max_degree = max_degree
  #         )
  #
  #
  #         ## if data-adaptive delta is true update delta for Q
  #         delta <- ind_gn_exp_estim$delta
  #
  #
  #         ## covariates are now exposures and baseline for Q
  #         covars <- c(a_names, w_names)
  #
  #         ## train Q with and without shifts
  #
  #         ind_qn_estim <- indiv_stoch_shift_est_Q(
  #           exposure = exposure,
  #           delta = delta,
  #           mu_learner = mu_learner,
  #           covars = covars,
  #           av = subset_av,
  #           at = subset_at,
  #           lower_bound = lower_bound,
  #           upper_bound = upper_bound,
  #           outcome_type = outcome_type
  #         )
  #
  #
  #         ## extract clever covariate
  #         Hn <- ind_gn_exp_estim$Hn_av
  #
  #
  #         ## do the TMLE update
  #         tmle_fit <- tmle_exposhift(
  #           data_internal = subset_av,
  #           delta = delta,
  #           Qn_scaled = ind_qn_estim$q_av,
  #           Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
  #           Hn = Hn,
  #           fluctuation = fluctuation,
  #           y = subset_av$y
  #         )
  #
  #         tmle_fit$call <- call
  #
  #         ## calculate the subpopulation intervention effect for the fold
  #         subpopulation_rank_shift_in_fold <- calc_final_ind_shift_param(
  #           tmle_fit,
  #           exposure,
  #           fold_k
  #         )
  #
  #         ## add in the delta shifted and the effect modifier for organization
  #         subpopulation_rank_shift_in_fold$Delta <- delta
  #         subpopulation_rank_shift_in_fold$Subgroup <- effect_mod_rule
  #
  #         ## save results as a list with fold, rank, exposure, modifier and
  #         ## level so we can extract results for pooling later
  #
  #         effect_mod_fold_results[[
  #           paste("Fold", ":", fold_k, "| Rank", ":", rank, "| Exposure", ":", exposure, "| Modifier", ":", effect_modifier, "| Level", ":", i)
  #         ]] <- list(
  #           "data" = subset_av,
  #           "Qn_scaled" = ind_qn_estim$q_av,
  #           "Hn" = Hn,
  #           "k_fold_result" = subpopulation_rank_shift_in_fold,
  #           "Delta" = delta,
  #           "Effect Mod Rule" = effect_mod_rule,
  #           "Exposure" = exposure
  #         )
  #       }
  #     }
  #
  #
  #     results_list <- list(
  #       effect_mod_fold_results
  #     )
  #
  #     names(results_list) <- c(
  #       "effect_mod_results"
  #     )
  #
  #     results_list
  #   },
  #   .options = furrr::furrr_options(seed = seed, packages = "EffectXshift")
  # )

  ## extract the results acros the folds
  k_fold_results <- do.call(rbind, purrr::map(fold_basis_results, c("k_fold_results")))
  exposure_shift_g <- unlist(purrr::map(fold_basis_results, c("exposure_shift_g")), recursive = FALSE)
  exposure_shift_q <- unlist(purrr::map(fold_basis_results, c("exposure_shift_q")), recursive = FALSE)
  g_estimates_region_v <- purrr::map(fold_basis_results, c("g_estimates_region_v"))
  g_estimates_region_vc <- purrr::map(fold_basis_results, c("g_estimates_region_vc"))
  q_estimates_region_v <- purrr::map(fold_basis_results, c("q_estimates_region_v"))
  q_estimates_region_vc <- purrr::map(fold_basis_results, c("q_estimates_region_vc"))
  data_region_v <- purrr::map(fold_basis_results, c("data_region_v"))
  data_region_vc <- purrr::map(fold_basis_results, c("data_region_vc"))
  data <- purrr::map(fold_basis_results, c("data"))



  ## calculate variable importance for each exposure
  pooled_exposure_results_list <- list()
  for (exposure in a_names) {
    q_fold_data_exposure <- exposure_shift_q[names(exposure_shift_q) == exposure]
    q_fold_data_exposure_combined <- do.call(rbind, q_fold_data_exposure)
    data_combined <- do.call(rbind, data)

    g_fold_data_exposure <- exposure_shift_g[names(exposure_shift_q) == exposure]
    g_fold_data_exposure_combined <- do.call(rbind, g_fold_data_exposure)

    tmle_fit <- tmle_exposhift(
      data_internal = data_combined,
      Qn_scaled = q_fold_data_exposure_combined,
      Hn = g_fold_data_exposure_combined,
      fluctuation = "standard",
      y = data_combined$y,
      delta = deltas[[exposure]]
    )

    indiv_shift_in_fold <- calc_final_ind_shift_param(
      tmle_fit = tmle_fit,
      exposure = exposure,
      fold_k = "Pooled TMLE"
    )

    pooled_exposure_results_list[[exposure]] <- indiv_shift_in_fold
  }

  pooled_exposure_results_df <- do.call(rbind, pooled_exposure_results_list)

  # get estimates for the v region:
  g_fold_data_region_combined <- do.call(rbind, g_estimates_region_v)
  q_fold_data_region_combined <- do.call(rbind, q_estimates_region_v)
  data_region_v_combined <- do.call(rbind, data_region_v)

  g_fold_data_region_combined <- g_fold_data_region_combined %>%
    mutate(across(everything(), ~ ifelse(. > hn_trunc_thresh, hn_trunc_thresh, .)))

  tmle_fit <- tmle_exposhift(
    data_internal = data_region_v_combined,
    Qn_scaled = q_fold_data_region_combined,
    Hn = g_fold_data_region_combined,
    fluctuation = "standard",
    y = data_region_v_combined$y,
    delta = mean(unlist(deltas))
  )

  indiv_shift_in_fold_v <- calc_final_ind_shift_param(
    tmle_fit = tmle_fit,
    exposure = "v",
    fold_k = "Pooled TMLE"
  )

  # get estimates for the v cregion:
  g_fold_data_region_combined <- do.call(rbind, g_estimates_region_vc)
  q_fold_data_region_combined <- do.call(rbind, q_estimates_region_vc)
  data_region_vc_combined <- do.call(rbind, data_region_vc)

  tmle_fit <- tmle_exposhift(
    data_internal = data_region_vc_combined,
    Qn_scaled = q_fold_data_region_combined,
    Hn = g_fold_data_region_combined,
    fluctuation = "standard",
    y = data_region_vc_combined$y,
    delta = mean(unlist(deltas))
  )

  indiv_shift_in_fold_vc <- calc_final_ind_shift_param(
    tmle_fit = tmle_fit,
    exposure = "vc",
    fold_k = "Pooled TMLE"
  )


  results_list <- list(
    "Effect Modification K-Fold Results" = k_fold_results,
    "Effect Modification Region V Pooled Results" = indiv_shift_in_fold_v,
    "Effect Modification Region V^c Pooled Results" = indiv_shift_in_fold_vc,
    "Marginal Shift Results" = pooled_exposure_results_df
  )

  return(results_list)
}
