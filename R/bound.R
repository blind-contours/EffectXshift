#' Bound Precision
#'
#' @details Bound values in the unit interval to machine precision in order to
#'  avoid numerical instability issues in downstream computation.
#'
#' @param vals \code{numeric} vector of values between 0 and 1 to be
#'  bounded within arbitrary machine precision. The most common use of this
#'  functionality is to avoid indeterminate or non-finite values after the
#'  application \code{stats::qlogis}.
#'
#' @importFrom assertthat assert_that
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#'  @export
bound_precision <- function(vals) {
  assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  vals[vals == 0] <- .Machine$double.neg.eps
  vals[vals == 1] <- 1 - .Machine$double.neg.eps
  return(vals)
}

###############################################################################

#' Bound Generalized Propensity Score
#'
#' @details Bound estimated values of the generalized propensity score (a
#'  conditional density) to avoid numerical instability issues arising from
#'  practical violations of the assumption of positivity.
#'
#' @param vals \code{numeric} vector of propensity score estimate values. Note
#'  that, for this parameter, the propensity score is (conditional) density and
#'  so it ought not be bounded from above.
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where the
#'  returned values are bounded such that the minimum is no lower than 1/n, for
#'  the sample size n.
#'  @export
bound_propensity <- function(vals) {
  # bound likelihood component g(a|w) away from 0 only
  propensity_bound <- 1 / length(vals)
  vals[vals < propensity_bound] <- propensity_bound
  return(vals)
}

###############################################################################

#' Transform values by scaling to the unit interval
#'
#' @details A transformation that scales an arbitrary set of input values to
#'  the unit interval. See \code{\link{scale_to_original}} for a corresponding
#'  backtransformation.
#'
#' @param vals A \code{numeric} vector corresponding to the observed values of
#'  the variable of interest, to be re-scaled to between 0 and 1.
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where the
#'  values are re-scaled to lie between 0 and 1.
#' @export
scale_to_unit <- function(vals) {
  # compute re-scaled value in interval [0,1]
  scaled_vals <- (vals - min(vals)) / (max(vals) - min(vals))
  return(scaled_vals)
}

###############################################################################

#' Transform values from the unit interval back to their original scale
#'
#' @details A back-transformation that returns values computed in the unit
#'  interval to their original scale. This is used in re-scaling updated TML
#'  estimates back to their natural scale. Undoes \code{\link{scale_to_unit}}.
#'
#' @param scaled_vals A \code{numeric} vector corresponding to re-scaled values
#'  in the unit interval, to be re-scaled to the original interval.
#' @param max_orig A \code{numeric} scalar value giving the maximum of the
#'  values on the original scale.
#' @param min_orig A \code{numeric} scalar value giving the minimum of the
#'  values on the original scale.
#'
#' @return A \code{numeric} vector of the same length as \code{scaled_vals},
#'  where the values are re-scaled to lie in their original/natural interval.
#' @export
scale_to_original <- function(scaled_vals, max_orig, min_orig) {
  scaled_orig <- scaled_vals * (max_orig - min_orig) + min_orig
  return(scaled_orig)
}
