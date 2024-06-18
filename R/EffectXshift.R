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
                         density_classification = FALSE,
                         max_depth = 1) {
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
        density_classification = density_classification,
        max_depth = max_depth
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
  g_fold_data_region_combined <- rbindlist(g_estimates_region_v)
  q_fold_data_region_combined <- rbindlist( q_estimates_region_v)
  data_region_v_combined <- rbindlist( data_region_v)

  tmle_fit <- tmle_exposhift(
    data_internal = data_region_v_combined,
    Qn_scaled = q_fold_data_region_combined,
    Qn_unscaled = scale_to_original(q_fold_data_region_combined, min_orig = min(data_region_v_combined$y), max_orig = max(data_region_v_combined$y)),
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
