#' Compute the Pooled Interaction Shift Parameter Estimate From the Fold
#' Specific Results
#'
#' @details Estimate the value of the pooled causal parameter alongside statistical
#'  inference for the parameter estimate based on the nuisance parameters from
#'  the fold specific results for the interaction parameter.
#'
#' @param intxn_shift_results A list of interaction results found across
#' the parallelized CV fold procedure
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#' @param a_names List of exposure names
#' @param w_names List of covariate names
#' @param z_names List of mediator names
#' @param fluctuation Type of fluctuation to be used
#'
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_intxn_shifts <- function(intxn_shift_results,
                                     estimator = c("tmle", "onestep"),
                                     fluc_mod_out = NULL,
                                     a_names,
                                     w_names,
                                     z_names,
                                     fluctuation,
                                     n_folds,
                                     rank = TRUE) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  names <- names(intxn_shift_results)

  if (rank == TRUE) {
    names <- gsub("^(Rank [0-9]+) :.*", "\\1", names)
  }else{
    names <- gsub("^.+?:", "", names)
  }
  k_fold_results_list <- list()
  pooled_results_list <- list()


  for (var_set in unique(names)) {
    var_set_results <- intxn_shift_results[stringr::str_detect(
      names(intxn_shift_results), var_set
    )]

    interaction_sets <- gsub("^.+?:", "", names(var_set_results))

    test <- unlist(var_set_results, recursive = FALSE)

    # Get clever covariate for each shift for each fold ----

    Hn <- test[stringr::str_detect(names(test), "Hn")]
    Hn_unlist <- unlist(Hn, recursive = FALSE)

    Hn_1 <- do.call(rbind, Hn_unlist[stringr::str_detect(
      names(Hn_unlist), "Hn1"
    )])
    Hn_2 <- do.call(rbind, Hn_unlist[stringr::str_detect(
      names(Hn_unlist), "Hn2"
    )])
    Hn_3 <- do.call(rbind, Hn_unlist[stringr::str_detect(
      names(Hn_unlist), "Hn3"
    )])

    # Get scaled Qn for each shift for each fold ----

    Qn <- test[stringr::str_detect(names(test), "Qn_scaled")]
    Qn_unlist <- unlist(Qn, recursive = FALSE)

    Qn_scaled_1 <- do.call(rbind, Qn_unlist[stringr::str_detect(
      names(Qn_unlist), "Qn_scaled1"
    )])
    Qn_scaled_2 <- do.call(rbind, Qn_unlist[stringr::str_detect(
      names(Qn_unlist), "Qn_scaled2"
    )])
    Qn_scaled_3 <- do.call(rbind, Qn_unlist[stringr::str_detect(
      names(Qn_unlist), "Qn_scaled3"
    )])

    data <- do.call(rbind, test[stringr::str_detect(
      names(test), "data"
    )])

    k_fold_results <- do.call(rbind, test[stringr::str_detect(
      names(test), "k_fold"
    )])


    deltas <- do.call(rbind, test[stringr::str_detect(
      names(test), "deltas"
    )])


    Qn_scaled <- list(Qn_scaled_1, Qn_scaled_2, Qn_scaled_3)
    Hn <- list(Hn_1, Hn_2, Hn_3)

    intxn_results_list <- list()


    for (i in 1:length(Hn)) {
      hn_estim <- Hn[[i]]
      qn_estim_scaled <- Qn_scaled[[i]]
      delta <- deltas[[i]]

      tmle_fit <- tmle_exposhift(
        data_internal = data,
        delta = delta,
        Qn_scaled = qn_estim_scaled,
        Hn = hn_estim,
        fluctuation = fluctuation,
        y = data$y,
      )

      intxn_results_list[[i]] <- tmle_fit
    }

    intxn_pooled <- calc_final_joint_shift_param(
      joint_shift_fold_results = intxn_results_list,
      rank = var_set,
      fold_k = "Pooled TMLE",
      deltas_updated = deltas,
      exposures = c("Var 1","Var 2","Joint","Interaction")
    )

    rownames(k_fold_results) <- NULL


    pooled_results_list[[var_set]] <- intxn_pooled
    k_fold_results_list[[var_set]] <- k_fold_results
  }


  return(list("k_fold_results" = k_fold_results_list, "pooled_results" = pooled_results_list))
}
