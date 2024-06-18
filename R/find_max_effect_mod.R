#' @title Identification of Maximal Effect Modifiers Using Stochastic Shift Interventions
#'
#' @description The `find_max_effect_mods` function identifies subpopulations with the maximum differential impact of stochastic shift interventions on exposures within a mixture. The method estimates the individual effects of shifting each exposure while controlling for other exposures and covariates, using ensemble machine learning and targeted maximum likelihood estimation (TMLE).
#' Once the individual intervention effects and influence curve estimates given an intervention on each exposure are derived, we then use a simple t-statistic partitioning algorithm to find the region with the maximum significant difference in intervention effects.
#'
#' @param data A \code{data.frame} containing all the variables needed for the analysis, including baseline covariates, exposures, and the outcome.
#' @param deltas A named \code{list} or \code{vector} specifying the shift in exposures to define the target parameter. Each element should correspond to an exposure variable specified in \code{a_names}, detailing the amount by which that exposure is to be shifted.
#' @param a_names A \code{character} vector specifying the names of the exposure variables within \code{data}.
#' @param w_names A \code{character} vector specifying the names of the covariate variables within \code{data}.
#' @param outcome The name of the outcome variable in \code{data}.
#' @param outcome_type A \code{character} string indicating the type of the outcome variable; either "continuous", "binary", or "count".
#' @param mu_learner A list of \code{\link[sl3]{Lrnr_sl}} learners specifying the ensemble machine learning models to be used for outcome prediction within the Super Learner framework.
#' @param g_learner A list of \code{\link[sl3]{Lrnr_sl}} learners specifying the ensemble machine learning models to be used for estimating the conditional density of the exposures.
#' @param top_n An \code{integer} specifying the number of top positive and negative effects to return.
#' @param seed An \code{integer} value to set the seed for reproducibility.
#' @param min_obs Minimum number of observations in a region to warrant a split.
#'
#' @return A list containing the top effects and interactions identified and estimated by the function. It includes elements for top positive and negative individual effects as well as top synergistic and antagonistic interactions.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(matrix(rnorm(100 * 10), ncol = 10))
#' names(data) <- c(paste0("X", 1:8), "exposure", "outcome")
#' deltas <- list(exposure = 0.1)
#' a_names <- "exposure"
#' w_names <- paste0("X", 1:8)
#' outcome <- "outcome"
#' outcome_type <- "continuous"
#' mu_learner <- list(sl3::Lrnr_mean$new(), sl3::Lrnr_glm$new())
#' top_n <- 3
#' seed <- 123
#'
#' results <- find_max_effect_mods(data, deltas, a_names, w_names, outcome, outcome_type, mu_learner, g_learner, top_n, seed, min_obs = 10)
#' print(results)
#' }
#'
#' @importFrom sl3 make_sl3_Task Lrnr_sl
#' @import dplyr
#' @import ranger
#' @importFrom data.table as.data.table setnames
#' @export

find_max_effect_mods <- function(at, av, deltas, a_names, w_names, outcome, outcome_type, mu_learner, g_learner, top_n = 3, seed, min_obs, fold, density_classification = TRUE, max_depth) {
  future::plan(future::sequential, gc = TRUE)
  set.seed(seed)

  # we need to impute any covariates before beginning the search algorithm
  at <- at %>%
    mutate(across(all_of(w_names), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  # we need to impute any covariates before beginning the search algorithm
  av <- av %>%
    mutate(across(all_of(w_names), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  # Initialize data frames for storing individual effects
  at_individual_effects_df <- list()
  at_influence_curve_df <- list()
  at_q_estimates <- list()
  at_g_estimates <- list()

  av_individual_effects_df <- list()
  av_influence_curve_df <- list()
  av_q_estimates <- list()
  av_g_estimates <- list()

  create_augmented_data <- function(at, av,  delta, var, covars) {
    n <- nrow(at)
    at_dup <- at[rep(1:n, each = 2), ]
    at_dup$intervention <- rep(c(0, 1), times = n)
    at_dup[[var]] <- ifelse(at_dup$intervention == 1, at_dup[[var]] + delta, at_dup[[var]])

    n <- nrow(av)
    av_dup <- av[rep(1:n, each = 2), ]
    av_dup$intervention <- rep(c(0, 1), times = n)
    av_dup[[var]] <- ifelse(av_dup$intervention == 1, av_dup[[var]] + delta, av_dup[[var]])
    return(list("av_dup" = av_dup, "at_dup" = at_dup))
  }

  # Function to estimate density ratio via classification
  estimate_density_ratio <- function(at, av, delta, var, covars, classifier) {
    augmented_data <- create_augmented_data(at, av, delta, var, covars)

    at_sl_task <- sl3::sl3_Task$new(
      data = augmented_data$at_dup,
      outcome = "intervention",
      covariates = covars,
      outcome_type = "binary"
    )

    av_sl_task <- sl3::sl3_Task$new(
      data = augmented_data$av_dup,
      outcome = "intervention",
      covariates = covars,
      outcome_type = "binary"
    )

    sl <- sl3::Lrnr_sl$new(
      learners = mu_learner,
      metalearner = sl3::Lrnr_nnls$new()
    )


    at_class_model <- suppressWarnings(suppressMessages(sl$train(at_sl_task)))

    # at predictions -----------
    at_class_model_preds <- at_class_model$predict(at_sl_task)
    av_class_model_preds <- at_class_model$predict(av_sl_task)

    # Compute the density ratios for shifted exposures
    at_u_t_shift <- at_class_model_preds[augmented_data$at_dup$intervention == 1]
    at_density_ratio_shift <- at_u_t_shift / (1 - at_u_t_shift)

    av_u_t_shift <- av_class_model_preds[augmented_data$av_dup$intervention == 1]
    av_density_ratio_shift <- av_u_t_shift / (1 - av_u_t_shift)

    # Compute the density ratios for unshifted exposures
    at_u_t_unshift <- at_class_model_preds[augmented_data$at_dup$intervention == 0]
    at_density_ratio_unshift <- 1

    av_u_t_unshift <- av_class_model_preds[augmented_data$av_dup$intervention == 0]
    av_density_ratio_unshift <- 1

    # Combine the density ratios into a data.table
    at_density_ratios <- data.table::as.data.table(
      cbind(at_density_ratio_unshift, at_density_ratio_shift)
    )

    data.table::setnames(at_density_ratios, c("noshift", "shift"))


    av_density_ratios <- data.table::as.data.table(
      cbind(av_density_ratio_unshift, av_density_ratio_shift)
    )

    data.table::setnames(av_density_ratios, c("noshift", "shift"))

    return(list(Hn_at = at_density_ratios, Hn_av = av_density_ratios))
  }



  # Calculate individual effects
  for (var in a_names) {

    ## get lower and upper bound so we don't shift too far
    lower_bound <- min(min(at[[var]]), min(at[[var]]))
    upper_bound <- max(max(at[[var]]), max(at[[var]]))

    if (density_classification) {
      ind_gn_exp_estim <- estimate_density_ratio(at = at, av = av, delta =  deltas[[var]], var = var, covars = c(var, w_names), classifier = mu_learner)
    } else {

    ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
      exposure = var,
      delta =  deltas[[var]],
      g_learner = g_learner,
      covars = w_names,
      av = av,
      at = at,
      adaptive_delta = FALSE,
      hn_trunc_thresh = hn_trunc_thresh,
      use_multinomial = FALSE,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      outcome_type = "continuous",
      density_type = "sl",
      n_bins = 4,
      max_degree = 10
    )
    }

    ## covariates are now exposures and baseline for Q
    covars <- c(a_names, w_names)

    ## train Q with and without shifts

    ind_qn_estim <- indiv_stoch_shift_est_Q(
      exposure = var,
      delta = deltas[[var]],
      mu_learner = mu_learner,
      covars = covars,
      av = av,
      at = at,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      outcome_type = outcome_type
    )


    ## extract clever covariate
    Hn_av <- ind_gn_exp_estim$Hn_av
    Hn_at <- ind_gn_exp_estim$Hn_at


    ## do the TMLE update for at
    at_tmle_fit <- tmle_exposhift(
      data_internal = at,
      delta =  deltas[[var]],
      Qn_scaled = ind_qn_estim$q_at,
      Qn_unscaled = scale_to_original(ind_qn_estim$q_at, min_orig = min(at$y), max_orig = max(at$y)),
      Hn = Hn_at,
      fluctuation = "standard",
      y = at$y
    )

    ## do the TMLE update for av
    av_tmle_fit <- tmle_exposhift(
      data_internal = av,
      delta =  deltas[[var]],
      Qn_scaled = ind_qn_estim$q_av,
      Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
      Hn = Hn_av,
      fluctuation = "standard",
      y = av$y
    )

    at_effect <- at_tmle_fit$qn_shift_star - at_tmle_fit$qn_noshift_star
    av_effect <- av_tmle_fit$qn_shift_star - av_tmle_fit$qn_noshift_star

    # Assuming you want to append the effect at the end of the dataframe.
    at_individual_effects_df[[var]] <- at_effect
    at_influence_curve_df[[var]] <- at_tmle_fit$eif
    at_q_estimates[[var]] <- ind_qn_estim$q_at
    at_g_estimates[[var]] <- ind_gn_exp_estim$Hn_at

    # Assuming you want to append the effect at the end of the dataframe.
    av_individual_effects_df[[var]] <- av_effect
    av_influence_curve_df[[var]] <- av_tmle_fit$eif
    av_q_estimates[[var]] <- ind_qn_estim$q_av
    av_g_estimates[[var]] <- ind_gn_exp_estim$Hn_av

  }

  at_individual_effects_df <- do.call(cbind, at_individual_effects_df)
  at_influence_curve_df <- do.call(cbind, at_influence_curve_df)

  calculate_p_value <- function(left_effects, left_ic, right_effects, right_ic) {
    left_mean <- mean(left_effects, na.rm = TRUE)
    right_mean <- mean(right_effects, na.rm = TRUE)
    left_se <- sqrt(var(left_ic, na.rm = TRUE) / length(left_ic))
    right_se <- sqrt(var(right_ic, na.rm = TRUE) / length(right_ic))
    t_stat <- (left_mean - right_mean) / sqrt(left_se^2 + right_se^2)
    p_value <- 2 * pt(-abs(t_stat), df = min(length(left_ic), length(right_ic)) - 1)
    return(p_value)
  }


  calculate_average <- function(region, w_names, var, outcome) {
    if (nrow(region) == 0) { # No data in this region
      return(list(average = Inf, n = 0)) # Use Inf as the default for average.
    } else {
      region_data <- region[, c(w_names, var, outcome)]
      # rf <- ranger(as.formula(paste(outcome, "~.")), data = region_data)

      region_average <- mean(region_data[[outcome]])
      return(list(average = region_average, n = nrow(region)))
    }
  }

  split_data <- function(data, split_var, split_point) {
    unique_values <- sort(unique(data[[split_var]]))
    # Check if the variable is binary (i.e., only contains 0 and 1)
    is_binary <- all(sort(unique_values) == c(0, 1))


    # Apply different ifelse conditions based on whether the variable is binary or not
    if (is_binary) {
      # If the variable is binary, apply the desired ifelse statement for binary variables
      left <- subset(data, data[[split_var]] == 0)
      right <- subset(data, data[[split_var]] == 1)
    } else {
      left <- subset(data, data[[split_var]] <= split_point)
      right <- subset(data, data[[split_var]] > split_point)
    }

    return(list(left = left, right = right))
  }

  clean_rule <- function(rule) {
    rule <- gsub("&\\s*$", "", rule)
    return(rule)
  }

  rules_to_dataframe <- function(rules_list) {
    rules_df <- data.frame(
      Average = numeric(),
      ParentAverage = numeric(),
      N = integer(),
      Rule = character(),
      Depth = integer(),
      stringsAsFactors = FALSE
    )

    for (rule in rules_list) {
      rules_df <- rbind(rules_df, data.frame(
        Average = rule$RegionMean,
        ParentAverage = rule$ParentMean,
        N = rule$N,
        Rule = rule$Rule,
        Depth = rule$Depth,
        stringsAsFactors = FALSE
      ))
    }

    return(rules_df)
  }

  find_best_split <- function(data, split_variables,
                              outcome, min_obs = 10, parent_p,
                              min_max) {
    best_split <- NULL
    min_p_value <- 1


    for (var in split_variables) {
      unique_values <- sort(unique(data[[var]]))
      # Check if the variable is binary (i.e., only contains 0 and 1)
      is_binary <- length(unique_values) == 2 && all(sort(unique_values) == c(0, 1))



      for (split_point in unique_values) {
        # Apply different ifelse conditions based on whether the variable is binary or not
        if (is_binary) {
          # If the variable is binary, apply the desired ifelse statement for binary variables
          data$split_indicator <- ifelse(data[[var]] == 1, 1, 0) # This is just an example, modify as needed
        } else {
          # If the variable is not binary, apply the original ifelse statement
          data$split_indicator <- ifelse(data[[var]] <= split_point, 1, 0)
        }

        # Check for minimum observations in each split
        if (sum(data$split_indicator) < min_obs || sum(!data$split_indicator)
        < min_obs) {
          next
        }

        left_data <- data[data$split_indicator == 1, ]
        right_data <- data[data$split_indicator == 0, ]

        left_effects <- left_data[[outcome]]
        right_effects <- right_data[[outcome]]

        left_ic <- left_data$ic
        right_ic <- right_data$ic

        p_value <- calculate_p_value(left_effects, left_ic, right_effects, right_ic)


        if (p_value < min_p_value & p_value <= parent_p) {
          min_p_value <- p_value
          best_split <- list(
            variable = var,
            point = split_point,
            p_value = min_p_value
          )
        }
      }
    }

    data$split_indicator <- NULL

    if (!is.null(best_split)) {
      return(best_split)
    } else {
      return(NULL) # No split found that reduces the average below the parent's
    }
  }

  recursive_split_all_rules <- function(data,
                                        split_variables,
                                        depth = 0,
                                        max_depth = 2,
                                        outcome,
                                        path = "",
                                        min_obs = min_obs,
                                        parent_p = Inf) {


    if (depth == max_depth || nrow(data) == 0) {
      current_rule <- list(
        "RegionMean" = mean(data[[outcome]], na.rm = TRUE), # Updated to use parent_stats
        "ParentMean" = parent_p,
        "N" = dim(data)[1], # Updated to use parent_stats
        "Rule" = clean_rule(path),
        "Depth" = depth
      )
      return(list(current_rule))
    }

    # Find the best split using the parent mean calculated earlier
    best_split <- find_best_split(data = data, split_variables = split_variables, outcome,
      min_obs = min_obs,
      parent_p = parent_p,
      min_max = min_max
    )

    if (is.null(best_split)) { # If no best split found, return the current path
      current_rule <- list(
        "RegionMean" = mean(data[[outcome]], na.rm = TRUE), # Updated to use parent_stats
        "ParentP" = parent_p,
        "N" = dim(data)[1], # Updated to use parent_stats
        "Rule" = clean_rule(path),
        "Depth" = depth
      )
      return(list(current_rule))
    }

    # Perform the split
    splits <- split_data(data = data, split_var = best_split$variable, split_point = best_split$point)

    unique_values <- sort(unique(data[[best_split$variable]]))
    # Check if the variable is binary (i.e., only contains 0 and 1)
    is_binary <- all(sort(unique_values) == c(0, 1))

    if (is_binary) {
      # Construct the rules for the left and right branches
      left_rule <- paste0(
        path, best_split$variable, " == ",
        0, " & "
      )
      right_rule <- paste0(
        path, best_split$variable, " == ",
        1, " & "
      )
    } else {
      left_rule <- paste0(
        path, best_split$variable, " <= ",
        best_split$point, " & "
      )
      right_rule <- paste0(
        path, best_split$variable, " > ",
        best_split$point, " & "
      )
    }

    left_rules <- recursive_split_all_rules(
      data = splits$left,
      split_variables = split_variables,
      depth = depth + 1,
      max_depth = max_depth,
      outcome,
      left_rule,
      min_obs,
      parent_p = best_split$p_value
    )
    right_rules <- recursive_split_all_rules(
      data = splits$right,
      split_variables = split_variables,
      depth = depth + 1,
      max_depth = max_depth,
      outcome,
      right_rule,
      min_obs,
      parent_p = best_split$p_value
    )

    # Combine the rules from the left and right branches and return them
    return(c(left_rules, right_rules))
  }

  # Go through the individual exposure effects and and find the region for each exposure that has the maximum difference

  effect_mod_results <- list()

  for (indiv_effect in 1:ncol(at_individual_effects_df)) {
    outcome <- at_individual_effects_df[, indiv_effect]
    ic <- at_influence_curve_df[, indiv_effect]
    mod_data <- as.data.frame(at)
    mod_data$y <- outcome
    mod_data$ic <- ic

    min_ave_tree_results <- recursive_split_all_rules(
      data = as.data.frame(mod_data),
      split_variables = w_names,
      max_depth = max_depth,
      outcome = "y",
      min_obs = min_obs
    )

    iteration_results <- as.data.frame(do.call(rbind, min_ave_tree_results))
    iteration_results$Exposure <- colnames(at_individual_effects_df)[indiv_effect]
    # iteration_results$exposure <- colnames(individual_effects_df)[indiv_effect]
    effect_mod_results[[colnames(at_individual_effects_df)[indiv_effect]]] <- iteration_results
  }

  # Calculate the difference in RegionMean values for each item and store them with their names
  region_mean_diffs <- sapply(effect_mod_results, function(x) abs(x$RegionMean[[1]] - x$RegionMean[[2]]))

  # Rank the items based on the differences
  ranked_names <- names(sort(region_mean_diffs, decreasing = TRUE))

  # Number of items you want to select

  # Select the top n items
  selected_items <- effect_mod_results[ranked_names[1:top_n]]
  results_list <- list()
  q_region_list <- list()
  g_region_list <- list()
  region_data_list <- list()

  for (item_index in 1:length(selected_items)) {
    item <- selected_items[[item_index]]
    item$Rank <- item_index
    selected_items[[item_index]] <- item
  }

  selected_items <- selected_items[[1]]

  for (row in 1:nrow(selected_items)) {
    item <- selected_items[row,]
    rule <- item$Rule
    exposure <- item$Exposure

    effect_indicator <- av %>%
      mutate(Indicator = eval(parse(text = paste0("(", rule, ")")))) %>%
      pull(Indicator)

    effect_indicator <- as.numeric(effect_indicator)

    effect_in_region <- mean(av_individual_effects_df[[exposure]][effect_indicator == 1], na.rm = TRUE)
    ic_in_region <- av_influence_curve_df[[exposure]][effect_indicator == 1]
    region_data <- av[effect_indicator == 1]

    av_q_estimates_region <- av_q_estimates[[exposure]][effect_indicator == 1]
    av_g_estimates_region <- av_g_estimates[[exposure]][effect_indicator == 1]

    # Compute standard error and confidence intervals
    variance_est <- var(ic_in_region, na.rm = TRUE) / length(ic_in_region)
    se_est <- sqrt(variance_est)
    CI <- calc_CIs(effect_in_region, se_est)

    Lower_CI <- CI[1]
    Upper_CI <- CI[2]

    # Store results in a data frame
    results <- data.frame(
      Exposure = exposure,
      Effect = effect_in_region,
      SE = se_est,
      `Lower CI` = Lower_CI,
      `Upper CI` = Upper_CI,
      Modifier = rule[[1]],
      Fold = fold
    )
    results_list[[row]] <- results
    g_region_list[[row]] <- av_g_estimates_region
    q_region_list[[row]] <- av_q_estimates_region
    region_data_list[[row]] <- region_data

  }

  results_df <- do.call(rbind, results_list)


  # Return the results
  return(list("K-fold_EM_results" = results_df,
              "av_q_estimates" = av_q_estimates,
              "av_hn_estimates" = av_g_estimates,
              "g_region_v" = g_region_list[[1]],
              "g_region_vc" = g_region_list[[2]],
              "q_region_v" = q_region_list[[1]],
              "q_region_vc" = q_region_list[[2]],
              "region_data_v" = region_data_list[[1]],
              "region_data_vc" = region_data_list[[2]],
              "data" = av))
}
