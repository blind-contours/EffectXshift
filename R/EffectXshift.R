#' @title EffectXshift
#'
#' @description Under a fixed shift to exposures identify using g-computation the joint shift of
#' pairwise exposures in a mixed exposure compared to the additive individual shifts. Positive values indicate
#' synergy and negative antagonism, get the top synergy and antagonism results and use CV-TMLE to efficiently
#' estimate the interaction target parameter.
#'
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#' baseline covariates. These variables are measured before exposures.
#' @param a \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple exposures.
#' @param z \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple mediators (optional).
#' @param y \code{numeric} vector of observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in exposures to
#' define the target parameter, with respect to the scale of the exposures (A). If adaptive_delta
#' is true, these values will be reduced.
#' @param var_sets A list specifying variable sets for deterministic EffectXshift usage.
#' Example: var_sets <- c("A_1", "A_1-Z_2") where the analyst provides variable sets
#' for exposures, exposure-mediator, or exposure-covariate relationships.
#' @param estimator The type of estimator to fit: \code{"tmle"} for targeted
#' maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param fluctuation Method used in the targeting step for TML estimation: "standard" or "weighted".
#' This determines where to place the auxiliary covariate in the logistic tilting regression.
#' @param pi_learner Learners for fitting Super Learner ensembles to densities via \pkg{sl3}.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param g_learner Learners for fitting Super Learner ensembles to the g-mechanism
#' g(A|W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param e_learner Learners for fitting Super Learner ensembles to the e-mechanism
#' g(A|Z,W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param zeta_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}..
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
#'
#' @return An S3 object of class \code{SuperNOVA} containing the results of the
#' procedure to compute a TML or one-step estimate of the counterfactual mean
#' under a modified treatment policy that shifts a continuous-valued exposure
#' by a scalar amount \code{delta}. These exposures are data-adaptively
#' identified using the CV-TMLE procedure.
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
                      var_sets = NULL,
                      pi_learner = NULL,
                      mu_learner = NULL,
                      g_learner = NULL,
                      e_learner = NULL,
                      zeta_learner = NULL,
                      n_folds = 2,
                      outcome_type = "continuous",
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 10,
                      adaptive_delta = FALSE,
                      discover_only = FALSE,
                      top_n = 1,
                      min_obs = 20) {
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

  if (is.null(pi_learner)) {
    sls <- create_sls()
    pi_learner <- sls$pi_learner
  }

  if (is.null(mu_learner)) {
    sls <- create_sls()
    mu_learner <- sls$mu_learner
  }

  if (is.null(zeta_learner)) {
    sls <- create_sls()
    zeta_learner <- sls$zeta_learner
  }

  if (is.null(g_learner)) {
    sls <- create_sls()
    g_learner <- sls$g_learner
  }

  if (is.null(e_learner)) {
    sls <- create_sls()
    e_learner <- sls$e_learner
  }

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

  data_internal <- data.table::data.table(w, a, y)
  `%notin%` <- Negate(`%in%`)

  if (outcome_type == "binary") {
    ## create the CV folds
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
  } else {
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
  }

  fold_basis_results <- furrr::future_map(unique(data_internal$folds),
                                          function(fold_k) {
                                            at <- data_internal[data_internal$folds != fold_k, ]
                                            av <- data_internal[data_internal$folds == fold_k, ]

                                            effect_mod_results <- find_max_effect_mods(
                                              data = at,
                                              deltas = deltas,
                                              a_names = a_names,
                                              w_names = w_names,
                                              outcome = "y",
                                              outcome_type = outcome_type,
                                              mu_learner = mu_learner,
                                              seed = seed,
                                              top_n = top_n,
                                              min_obs = min_obs
                                            )

                                          },
                                          .options = furrr::furrr_options(seed = seed, packages = "EffectXshift")
  )

  effect_mod_fold_results <- list()

  fold_effectXshift_results <- furrr::future_map(
    unique(data_internal$folds), function(fold_k) {

      fold_intxn_results <- fold_basis_results[[fold_k]]

      for (rank in 1:length(fold_intxn_results)) {
        rank_results <- fold_intxn_results[[rank]]


      for (i in 1:nrow(rank_results)) {

        rank_row <- rank_results[i,]

        exposure <- rank_row$Exposure
        effect_mod_rule <- unlist(rank_row$Rule)
        rank <- rank_row$Rank

        effect_modifier <- Filter(function(item) grepl(item, effect_mod_rule), w_names)


        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]

        delta <- deltas[[exposure]]

        lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
        upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))

        subset_at <- subset(at, eval(parse(text = effect_mod_rule)))
        subset_av <- subset(av, eval(parse(text = effect_mod_rule)))


        ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
          exposure = exposure,
          delta = delta,
          g_learner = pi_learner,
          covars = w_names,
          av = subset_av,
          at = subset_at,
          adaptive_delta = adaptive_delta,
          hn_trunc_thresh = hn_trunc_thresh,
          use_multinomial = FALSE,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = "continuous",
          density_type = "sl",
          n_bins = n_bins,
          max_degree = max_degree
        )

        delta <- ind_gn_exp_estim$delta

        covars <- c(a_names, w_names)

        ind_qn_estim <- indiv_stoch_shift_est_Q(
          exposure = exposure,
          delta = delta,
          mu_learner = mu_learner,
          covars = covars,
          av = subset_av,
          at = subset_at,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = outcome_type
        )

        Hn <- ind_gn_exp_estim$Hn_av

        tmle_fit <- tmle_exposhift(
          data_internal = subset_av,
          delta = delta,
          Qn_scaled = ind_qn_estim$q_av,
          Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
          Hn = Hn,
          fluctuation = fluctuation,
          y = subset_av$y
        )

        tmle_fit$call <- call

        subpopulation_rank_shift_in_fold <- calc_final_ind_shift_param(
          tmle_fit,
          exposure,
          fold_k
        )

        subpopulation_rank_shift_in_fold$Delta <- delta
        subpopulation_rank_shift_in_fold$Subgroup <- effect_mod_rule

        effect_mod_fold_results[[
          paste("Fold", ":", fold_k, "| Rank", ":", rank, "| Exposure", ":", exposure, "| Modifier", ":", effect_modifier, "| Level", ":", i )
        ]] <- list(
          "data" = subset_av,
          "Qn_scaled" = ind_qn_estim$q_av,
          "Hn" = Hn,
          "k_fold_result" = subpopulation_rank_shift_in_fold,
          "Delta" = delta,
          "Effect Mod Rule" = effect_mod_rule,
          "Exposure" = exposure
        )

      }
      }


      results_list <- list(
        effect_mod_fold_results
      )

      names(results_list) <- c(
        "effect_mod_results"

      )

      results_list
    },
    .options = furrr::furrr_options(seed = seed, packages = "EffectXshift")
  )

  ranked_effect_mod_results <- purrr::map(fold_effectXshift_results, c("effect_mod_results"))

  ranked_effect_mod_results <- unlist(ranked_effect_mod_results, recursive = FALSE)


  ## calculate pooled results based on rank
  effect_modification_results <- calc_pooled_indiv_shifts(
    indiv_shift_results = ranked_effect_mod_results,
    estimator = estimator,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = TRUE
  )


  results_list <- list(
    "Effect Modification K-Fold Results" = effect_modification_results$k_fold_results,
    "Effect Modification Pooled Results" = effect_modification_results$pooled_results
  )

  return(results_list)
}
