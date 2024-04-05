#' @title Data-adaptive Discovery of Interactions Based on Joint vs. Individual Shift Interventions
#'
#' @description The `find_max_effect_mods` function provides a g-computation approach to finding interactions.
#' This implementation fits an SL and then predicts two way joint shifts in exposures and compares this to individual shifts,
#' ranks the interactions based on the difference and then this becomes the interactions we want to estimate
#' using CV-TMLE.
#'
#' @param data A \code{data.frame} containing all the variables needed for the analysis, including
#' baseline covariates, exposures, and the outcome.
#' @param deltas A named \code{list} or \code{vector} specifying the shift in exposures to define the target parameter.
#' Each element should correspond to an exposure variable specified in \code{a_names}, detailing the amount
#' by which that exposure is to be shifted.
#' @param a_names A \code{character} vector specifying the names of the exposure variables within \code{data}.
#' @param w_names A \code{character} vector specifying the names of the covariate variables within \code{data}.
#' @param outcome The name of the outcome variable in \code{data}.
#' @param outcome_type A \code{character} string indicating the type of the outcome variable; either "continuous",
#' "binary", or "count".
#' @param mu_learner A list of \code{\link[sl3]{Lrnr_sl}} learners specifying the ensemble machine learning models to be used
#' for outcome prediction within the Super Learner framework.
#' @param top_n An \code{integer} specifying the number of top positive and negative effects to return.
#' @param seed An \code{integer} value to set the seed for reproducibility.
#'
#' @return A list containing the top effects and interactions identified and estimated by the function.
#' It includes elements for top positive and negative individual effects as well as top synergistic and antagonistic interactions.
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
#' results <- find_synergy_antagonism(data, deltas, a_names, w_names, outcome, outcome_type, mu_learner, top_n, seed)
#' print(results)
#' }
#'
#' @importFrom sl3 make_sl3_Task Lrnr_sl
#' @import dplyr
#' @import ranger
#' @export

find_max_effect_mods <- function(data, deltas, a_names, w_names, outcome, outcome_type, mu_learner, top_n = 3, seed, min_obs) {
  future::plan(future::sequential, gc = TRUE)
  set.seed(seed)

  # we need to impute any covariates before beginning the search algorithm
  data <- data %>%
    mutate(across(all_of(w_names), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  # Prepare Super Learner Task with specified data, exposures, and outcome
  task <- sl3::make_sl3_Task(
    data = data,
    covariates = c(a_names, w_names), # w_names assumed to be predefined
    outcome = outcome,
    outcome_type = outcome_type
  )

  # Train Super Learner
  sl <- sl3::Lrnr_sl$new(learners = mu_learner)
  sl_fit <- sl$train(task)

  # Initialize data frames for storing individual effects
  individual_effects_df <- list()

  # Calculate individual effects
  for (var in a_names) {
    shifted_data <- as.data.frame(data)
    shifted_data[[var]] <- shifted_data[[var]] + deltas[[var]]
    predictions <- sl_fit$predict(sl3::make_sl3_Task(data = shifted_data, covariates = c(a_names, w_names), outcome = outcome))
    effect <- predictions - data[[outcome]]

    # Assuming you want to append the effect at the end of the dataframe.
    individual_effects_df[[var]] <- effect
  }

  individual_effects_df <- do.call(cbind, individual_effects_df)


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

  find_best_split <- function(data, split_variables, w_names,
                              outcome, min_obs = 10, parent_average,
                              min_max) {
    best_split <- NULL
    min_average <- Inf
    max_average <- -Inf

    for (var in split_variables) {
      unique_values <- sort(unique(data[[var]]))
      # Check if the variable is binary (i.e., only contains 0 and 1)
      is_binary <- all(sort(unique_values) == c(0, 1))


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

        covars <- setdiff(w_names, var)

        # model <- ranger(
        #   formula = as.formula(paste(
        #     outcome, "~ split_indicator +",
        #     paste(covars, collapse = "+")
        #   )),
        #   data = data, num.trees = 400
        # )



        # Calculate the average predictions for left and right
        left_average <- mean(data[[outcome]][data$split_indicator == 1])
        right_average <- mean(data[[outcome]][data$split_indicator == 0])

        abs_ate <- abs(right_average - left_average)

        if (abs_ate > max_average && abs_ate > parent_average) {
          max_average <- abs_ate
          best_split <- list(
            variable = var,
            point = split_point,
            average = max_average
          )
        }
      }
    }

    # Clean up the temporary variable
    data$split_indicator <- NULL

    if (!is.null(best_split)) {
      return(best_split)
    } else {
      return(NULL) # No split found that reduces the average below the parent's
    }
  }

  recursive_split_all_rules <- function(data,
                                        split_variables,
                                        w_names,
                                        depth = 0,
                                        max_depth = 2,
                                        outcome,
                                        path = "",
                                        min_obs = min_obs,
                                        parent_mean = ifelse(min_max == "min", Inf, -Inf),
                                        min_max) {
    # Calculate the parent mean before attempting to find the best split
    region_stats <- calculate_average(
      region = data,
      w_names = w_names,
      var = split_variables,
      outcome = outcome
    )

    if (depth == max_depth || nrow(data) == 0) {
      current_rule <- list(
        "Value" = mean(data[[outcome]], na.rm = TRUE),
        "RegionMean" = region_stats$average, # Updated to use parent_stats
        "ParentMean" = parent_mean,
        "N" = region_stats$n, # Updated to use parent_stats
        "Rule" = clean_rule(path),
        "Depth" = depth
      )
      return(list(current_rule))
    }

    # Find the best split using the parent mean calculated earlier
    best_split <- find_best_split(data, split_variables, w_names, outcome,
      min_obs,
      parent_average = parent_mean,
      min_max = min_max
    )

    if (is.null(best_split)) { # If no best split found, return the current path
      current_rule <- list(
        "Value" = mean(data[[outcome]], na.rm = TRUE),
        "RegionMean" = region_stats$average, # Updated to use parent_stats
        "ParentMean" = parent_mean,
        "N" = region_stats$n, # Updated to use parent_stats
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
      w_names = w_names,
      depth = depth + 1,
      max_depth = max_depth,
      outcome,
      left_rule,
      min_obs,
      parent_mean = best_split$average,
      min_max = min_max
    )
    right_rules <- recursive_split_all_rules(
      data = splits$right,
      split_variables = split_variables,
      w_names = w_names,
      depth = depth + 1,
      max_depth = max_depth,
      outcome,
      right_rule,
      min_obs,
      parent_mean = best_split$average,
      min_max = min_max
    )

    # Combine the rules from the left and right branches and return them
    return(c(left_rules, right_rules))
  }

  # Go through the individual exposure effects and and find the region for each exposure that has the maximum difference

  effect_mod_results <- list()

  for (indiv_effect in 1:ncol(individual_effects_df)) {
    outcome <- individual_effects_df[, indiv_effect]
    mod_data <- as.data.frame(data)
    mod_data$y <- outcome
    exposures_to_adj <- a_names[-indiv_effect]

    min_ave_tree_results <- recursive_split_all_rules(
      data = as.data.frame(mod_data),
      split_variables = w_names,
      w_names = c(w_names, exposures_to_adj),
      max_depth = 1,
      outcome = "y",
      min_max = "max",
      min_obs = min_obs
    )

    iteration_results <- as.data.frame(do.call(rbind, min_ave_tree_results))
    iteration_results$Exposure <- colnames(individual_effects_df)[indiv_effect]
    # iteration_results$exposure <- colnames(individual_effects_df)[indiv_effect]
    effect_mod_results[[colnames(individual_effects_df)[indiv_effect]]] <- iteration_results
  }

  # Calculate the difference in RegionMean values for each item and store them with their names
  region_mean_diffs <- sapply(effect_mod_results, function(x) abs(x$RegionMean[[1]] - x$RegionMean[[2]]))

  # Rank the items based on the differences
  ranked_names <- names(sort(region_mean_diffs, decreasing = TRUE))

  # Number of items you want to select

  # Select the top n items
  selected_items <- effect_mod_results[ranked_names[1:top_n]]

  for (item_index in 1:length(selected_items)) {
    item <- selected_items[[item_index]]
    item$Rank <- item_index
    selected_items[[item_index]] <- item
  }


  # Return the results
  return(selected_items)
}
