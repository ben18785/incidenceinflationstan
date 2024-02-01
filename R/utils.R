#' Generate weights series
#'
#' @param t_max maximum time to generate weights until
#' @inheritParams gamma_discrete_pmf
#'
#' @return a vector of normalised weights of length t_max
weights_series <- function(t_max, serial_parameters) {
  day_series <- seq(1, t_max, 1)
  w <- purrr::map_dbl(day_series, ~gamma_discrete_pmf(., serial_parameters))
  # due to truncation, normalise series
  w <- w / sum(w)
  w
}

#' Gives cases expected given history of cases and Rt
#'
#' @inheritParams state_process_logp
#' @param weights a vector of normalised weights
#'
#' @return an expected number of cases
expected_cases <- function(Rt, weights, cases_history) {
  if(length(weights) != length(cases_history))
    stop("weights and history of cases must be same length.")
  Rt * sum(weights * cases_history)
}

#' Reparameterisation of gamma cdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#'
#' @return a cumulative density
pgamma_mean_sd <- function(x, mu, sigma) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::pgamma(x, shape, rate)
}

#' Reparameterisation of gamma pdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#' @param ... other arguments passed to method
#' @return a density
dgamma_mean_sd <- function(x, mu, sigma, ...) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::dgamma(x, shape, rate, ...)
}

#' Reparameterisation of gamma inverse-cdf using mean and sd
#'
#' @param x a value to be evaluated at
#' @param mu a mean
#' @param sigma a standard deviation
#'
#' @return a positive value
qgamma_mean_sd <- function(x, mu, sigma) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::qgamma(x, shape, rate)
}

#' Discrete gamma probability mass function
#'
#' @param day day to evaluate pmf
#' @param serial_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the serial interval distribution
#'
#' @return a probability
gamma_discrete_pmf <- function(day, serial_parameters){
  delay_mean <- serial_parameters$mean
  delay_sd <- serial_parameters$sd
  pgamma_mean_sd(day + 0.5, delay_mean, delay_sd) -
    pgamma_mean_sd(day - 0.5, delay_mean, delay_sd)
}

check_parameter_names <- function(reporting_parameters, serial_parameters) {
  correct_names_r <- c("location", "scale")
  correct_names_s <- c("mean", "sd")
  if(!all(correct_names_r %in% names(reporting_parameters)))
    stop("reporting_parameters must contain 'location' and 'scale'.")
  if(!all(correct_names_s %in% names(serial_parameters)))
    stop("serial_parameters must contain 'mean' and 'sd'.")
}
