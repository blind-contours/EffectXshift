#' Compute the Pooled Shift Parameter Estimate From the Fold Specific Results
#'
#' @details Estimate the value of the pooledcausal parameter alongside
#' statistical inference for the parameter estimate based on the nuisance
#' parameters from the fold specific results.
#'
#' @param indiv_shift_results List of individual shift results from the
#' parallelized CV estimation
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  estimator.
#' @param fluc_mod_out An object giving values of the logistic tilting model
#'  for targeted minimum loss estimation. This type of object should be the
#'  output of the internal routines to perform this step of the TML estimation
#'  procedure, as given by \code{\link{fit_fluctuation}}.
#' @param fluctuation Type of fluctuation to use
#' @param n_folds Number of folds used in the CV.
#' @param rank If we are pooling by rank
#'
#' @importFrom stats var
#'
#' @return A \code{list} containing the parameter estimate, estimated variance
#'  based on the efficient influence function (EIF), the estimate of the EIF
#'  incorporating inverse probability of censoring weights, and the estimate of
#'  the EIF without the application of such weights.
calc_pooled_indiv_shifts <- function(indiv_shift_results,
                                     estimator = c("tmle", "onestep"),
                                     fluc_mod_out = NULL,
                                     fluctuation,
                                     n_folds,
                                     rank = TRUE) {
  # set TMLE as default estimator type
  estimator <- match.arg(estimator)

  k_fold_results_list <- list()
  pooled_fold_results_list <- list()

  names <- names(indiv_shift_results)

  ## we want to extract data for the ranks and levels which we will pool on
  pattern <- "Fold : (\\d) \\| Rank : (\\d) \\| Exposure : ([^\\|]+) \\| Modifier : ([^\\|]+) \\| Level : (\\d)"
  ids <- sapply(names, function(name) {
    matches <- stringr::str_match(name, pattern)
    if (!is.na(matches[1, 2])) { # Check if there was a match
      paste(matches[, 3], matches[, 6], sep = "_")
    } else {
      NA
    }
  })

  ## for the rank level combinations:
  for (id in unique(ids)) {
    # Filter results for the current id
    id_results <- indiv_shift_results[ids == id]

    ## get the relevant data across the folds and create dfs for each

    Hn <- do.call(rbind, lapply(id_results, function(x) x$Hn))
    Qn_scaled <- do.call(rbind, lapply(id_results, function(x) x$Qn_scaled))
    data <- do.call(rbind, lapply(id_results, function(x) x$data))
    k_fold_results <- do.call(rbind, lapply(id_results, function(x) x$k_fold))
    deltas <- do.call(rbind, lapply(id_results, function(x) x$Delta))
    rule <- do.call(rbind, lapply(id_results, function(x) x$`Effect Mod Rule`))

    ## apply the TMLE update to these across all the folds together

    tmle_fit <- tmle_exposhift(
      delta = mean(deltas),
      data_internal = data,
      Qn_scaled = Qn_scaled,
      Hn = Hn,
      fluctuation = fluctuation,
      y = data$y
    )

    ## apply the function to put together the single shift results

    indiv_shift_in_fold <- calc_final_ind_shift_param(
      tmle_fit = tmle_fit,
      exposure = id,
      fold_k = "Pooled TMLE"
    )

    indiv_shift_in_fold$Delta <- mean(deltas)

    pooled_fold_results_list[[id]] <- indiv_shift_in_fold
  }


  for (fold in names(indiv_shift_results)) {
    id_results <- indiv_shift_results[fold]

    Hn <- do.call(rbind, lapply(id_results, function(x) x$Hn))
    Qn_scaled <- do.call(rbind, lapply(id_results, function(x) x$Qn_scaled))
    data <- do.call(rbind, lapply(id_results, function(x) x$data))
    k_fold_results <- do.call(rbind, lapply(id_results, function(x) x$k_fold))
    deltas <- do.call(rbind, lapply(id_results, function(x) x$Delta))
    rule <- do.call(rbind, lapply(id_results, function(x) x$`Effect Mod Rule`))


    tmle_fit <- tmle_exposhift(
      data_internal = data,
      Qn_scaled = Qn_scaled,
      Hn = Hn,
      fluctuation = fluctuation,
      y = data$y,
      delta = mean(deltas)
    )

    indiv_shift_in_fold <- calc_final_ind_shift_param(
      tmle_fit = tmle_fit,
      exposure = fold,
      fold_k = "K-fold TMLE"
    )

    indiv_shift_in_fold$Delta <- mean(deltas)
    indiv_shift_in_fold$Covariate_region <- rule

    k_fold_results_list[[fold]] <- indiv_shift_in_fold
  }

  fold_results_list <- list()

  for (fold in 1:n_folds) {
    # This pattern matches items for the current fold
    pattern <- sprintf("Fold : %d", fold)

    # Find keys (names) that match the current fold
    fold_keys <- names(k_fold_results_list)[grepl(pattern, names(k_fold_results_list))]

    # Extract the data frames for the current fold
    fold_data <- k_fold_results_list[fold_keys]

    # Combine all data frames for the current fold into a single data frame
    combined_fold_data <- do.call("rbind", fold_data)

    # Assign the combined data frame to the fold_results_list
    fold_results_list[[paste("Fold", fold)]] <- combined_fold_data
  }

  return(list("k_fold_results" = fold_results_list, "pooled_results" = pooled_fold_results_list))
}
