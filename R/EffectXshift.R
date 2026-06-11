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
#' @param estimator The type of estimator to fit: \code{"tmle"} (default) for
#' targeted maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param g_learner Learners for fitting the exposure mechanism g(A|W) via \pkg{sl3}.
#' Used to estimate the density ratio for continuous exposures when
#' \code{density_classification = FALSE}.
#' @param n_folds Number of folds to use in cross-validation, default is 2.
#' @param outcome_type Data type of the outcome, default is "continuous".
#' @param parallel Whether to parallelize across cores (default: TRUE).
#' @param parallel_type Type of parallelization to use if parallel is TRUE:
#' "multi_session" (default), "multicore", or "sequential".
#' @param num_cores Number of CPU cores to use in parallelization (default: 2).
#' @param seed \code{numeric} seed value to be passed to all functions.
#' @param hn_trunc_thresh Truncation level for the clever covariate (default: 10).
#' @param adaptive_delta If TRUE, reduces the user-specified delta until
#' the Hn calculated for a shift does not have any observation greater
#' than hn_trunc_thresh (default: FALSE).
#' @param min_obs Minimum number of observations allowed in a covariate region.
#' @param top_n top number of effect modifier-exposure pairs to estimate.
#' @param density_classification If TRUE, estimate the exposure density ratio
#' via a classification reparameterization rather than direct conditional
#' density estimation (default: FALSE). Continuous-exposure path only.
#' @param max_depth Maximum depth of the partitioning tree used to discover
#' effect-modification regions (default: 1).
#' @param rct If TRUE, run the randomized-trial workflow for a single binary
#' exposure (estimate the subject-level treatment effect and find the oracle
#' effect-modification region). If FALSE (default), run the mixed-exposure
#' stochastic-shift workflow.
#' @param rct_type Only used when \code{rct = TRUE}. Either \code{"ate"} (default),
#' which targets the subject-level ATE \eqn{Q(1, W) - Q(0, W)}, or \code{"incps"},
#' an incremental propensity-score shift from \eqn{\alpha} to \eqn{\alpha + \delta}.
#' @param pval_thresh p-value threshold for accepting a partition split when
#' discovering effect-modification regions in the RCT workflow (default: 0.1).
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
                         estimator = c("tmle", "onestep"),
                         fluctuation = c("standard", "weighted"),
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
                         max_depth = 1,
                         rct = FALSE,
                         rct_type = c("ate", "incps"),
                         pval_thresh = 0.1) {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  rct_type <- match.arg(rct_type)
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

      if (rct == FALSE) {
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
      } else{
        effect_mod_results <- find_max_effect_mods_rct(
          at = at,
          av = av,
          delta = as.numeric(deltas),
          a_name = a_names,
          w_names = w_names,
          outcome = "y",
          outcome_type = outcome_type,
          mu_learner = mu_learner,
          seed = seed,
          top_n = top_n,
          min_obs = min_obs,
          fold = fold_k,
          max_depth = max_depth,
          rct_type = rct_type,
          pval_thresh = pval_thresh
        )



      }



      k_fold_results <- effect_mod_results$K_fold_EM_results
      exposures_shift_q_results <- effect_mod_results$av_q_estimates
      exposures_shift_g_results <- effect_mod_results$av_hn_estimates
      g_region_v <- effect_mod_results$g_region_v
      g_region_vc <- effect_mod_results$g_region_vc
      q_region_v <- effect_mod_results$q_region_v
      q_region_vc <- effect_mod_results$q_region_vc
      data_region_v <- effect_mod_results$data_region_v
      data_region_vc <- effect_mod_results$data_region_vc
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


  ## extract the results across the folds
  k_fold_results <- do.call(rbind, purrr::map(fold_basis_results, c("k_fold_results")))

  # Combine the region-based data from each fold
  g_estimates_region_v  <- purrr::map(fold_basis_results, c("g_estimates_region_v"))
  g_estimates_region_vc <- purrr::map(fold_basis_results, c("g_estimates_region_vc"))
  q_estimates_region_v  <- purrr::map(fold_basis_results, c("q_estimates_region_v"))
  q_estimates_region_vc <- purrr::map(fold_basis_results, c("q_estimates_region_vc"))
  data_region_v  <- purrr::map(fold_basis_results, c("data_region_v"))
  data_region_vc <- purrr::map(fold_basis_results, c("data_region_vc"))
  data <- purrr::map(fold_basis_results, c("data"))

  ## Are we in a scenario with multiple exposures or a single binary exposure?
  if (length(a_names) > 1) {
    # -------------------------------------------------
    # CASE 1: Multiple exposures in a mixture scenario
    # -------------------------------------------------

    # 1) Combine "exposure_shift_q" and "exposure_shift_g"
    exposure_shift_q <- unlist(purrr::map(fold_basis_results, c("exposure_shift_q")), recursive = FALSE)
    exposure_shift_g <- unlist(purrr::map(fold_basis_results, c("exposure_shift_g")), recursive = FALSE)

    # 2) For each exposure, bind them across folds, then run tmle_exposhift
    pooled_exposure_results_list <- list()
    for (exposure in a_names) {
      # subset the relevant fold results for this exposure
      q_fold_data_exposure <- exposure_shift_q[names(exposure_shift_q) == exposure]
      q_fold_data_exposure_combined <- do.call(rbind, q_fold_data_exposure)

      g_fold_data_exposure <- exposure_shift_g[names(exposure_shift_q) == exposure]
      g_fold_data_exposure_combined <- do.call(rbind, g_fold_data_exposure)

      # Also combine the underlying data for these folds
      data_combined <- do.call(rbind, data)

      # run tmle_exposhift
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

    # Combine results
    pooled_exposure_results_df <- do.call(rbind, pooled_exposure_results_list)

    # Then do region V, region V^c for the aggregated data
    g_fold_data_region_combined <- data.table::rbindlist(g_estimates_region_v)
    q_fold_data_region_combined <- data.table::rbindlist(q_estimates_region_v)
    data_region_v_combined      <- data.table::rbindlist(data_region_v)

    tmle_fit_v <- tmle_exposhift(
      data_internal = data_region_v_combined,
      Qn_scaled = q_fold_data_region_combined,
      # optional if you have scale_to_original or other transformations
      Hn = g_fold_data_region_combined,
      fluctuation = "standard",
      y = data_region_v_combined$y,
      delta = mean(unlist(deltas))
    )
    indiv_shift_in_fold_v <- calc_final_ind_shift_param(
      tmle_fit = tmle_fit_v,
      exposure = "v",
      fold_k = "Pooled TMLE"
    )

    # complement
    g_fold_data_region_combined_vc <- data.table::rbindlist(g_estimates_region_vc)
    q_fold_data_region_combined_vc <- data.table::rbindlist(q_estimates_region_vc)
    data_region_vc_combined        <- data.table::rbindlist(data_region_vc)

    tmle_fit_vc <- tmle_exposhift(
      data_internal = data_region_vc_combined,
      Qn_scaled = q_fold_data_region_combined_vc,
      Hn = g_fold_data_region_combined_vc,
      fluctuation = "standard",
      y = data_region_vc_combined$y,
      delta = mean(unlist(deltas))
    )
    indiv_shift_in_fold_vc <- calc_final_ind_shift_param(
      tmle_fit = tmle_fit_vc,
      exposure = "vc",
      fold_k = "Pooled TMLE"
    )

    results_list <- list(
      "Effect Modification K-Fold Results"           = k_fold_results,
      "Effect Modification Region V Pooled Results"  = indiv_shift_in_fold_v,
      "Effect Modification Region V^c Pooled Results" = indiv_shift_in_fold_vc,
      "Marginal Shift Results"                       = pooled_exposure_results_df
    )

    return(results_list)

  } else {
    # -------------------------------------------------
    # CASE 2: Single binary exposure in an RCT
    # -------------------------------------------------
    # find_max_effect_mods_rct() targets a subject-level treatment effect
    # (the ATE Q(1, W) - Q(0, W), or an incremental propensity shift) on the
    # training folds, discovers the region V maximizing the differential effect,
    # then evaluates V and V^c on the held-out validation folds, returning the
    # subject-level effect and its influence-curve contribution per observation.
    #
    # Pooling these validation contributions across folds is the CV-TMLE
    # estimate of the oracle parameter:
    #   psi_V    = E[effect | W in V]
    #   psi_V^c  = E[effect | W in V^c]
    #   contrast = psi_V - psi_V^c   (the oracle effect-modification parameter)

    data_region_v_combined  <- data.table::rbindlist(data_region_v)
    data_region_vc_combined <- data.table::rbindlist(data_region_vc)

    # Pooled estimate, SE, and Wald CI for a region using the per-subject
    # influence curve (stored in the `ic` column by find_max_effect_mods_rct).
    pool_region <- function(region_data, region_label) {
      eff <- region_data$Effect
      ic  <- region_data$ic
      n   <- length(eff)
      psi <- mean(eff, na.rm = TRUE)
      var_est <- stats::var(ic, na.rm = TRUE) / n
      se  <- sqrt(var_est)
      ci  <- calc_CIs(psi, se)
      p_value <- 2 * stats::pnorm(abs(psi / se), lower.tail = FALSE)
      data.frame(
        Region     = region_label,
        Psi        = psi,
        SE         = se,
        `Lower CI` = ci[1],
        `Upper CI` = ci[2],
        `P-value`  = p_value,
        N          = n,
        check.names = FALSE
      )
    }

    region_v_res  <- pool_region(data_region_v_combined, "V")
    region_vc_res <- pool_region(data_region_vc_combined, "V^c")

    # Oracle contrast V vs V^c. The two regions are disjoint sets of subjects,
    # so the contrast variance is the sum of the region variances.
    contrast_psi <- region_v_res$Psi - region_vc_res$Psi
    contrast_se  <- sqrt(region_v_res$SE^2 + region_vc_res$SE^2)
    contrast_ci  <- calc_CIs(contrast_psi, contrast_se)
    contrast_p   <- 2 * stats::pnorm(abs(contrast_psi / contrast_se),
                                     lower.tail = FALSE)

    contrast_res <- data.frame(
      Region     = "V - V^c",
      Psi        = contrast_psi,
      SE         = contrast_se,
      `Lower CI` = contrast_ci[1],
      `Upper CI` = contrast_ci[2],
      `P-value`  = contrast_p,
      N          = region_v_res$N + region_vc_res$N,
      check.names = FALSE
    )

    pooled_results <- rbind(region_v_res, region_vc_res, contrast_res)
    rownames(pooled_results) <- NULL

    results_list <- list(
      "Effect Modification K-Fold Results" = k_fold_results,
      "Pooled Region Effects"              = pooled_results,
      "Region V Data"                      = data_region_v_combined,
      "Region V^c Data"                    = data_region_vc_combined
    )

    return(results_list)
  }
}

