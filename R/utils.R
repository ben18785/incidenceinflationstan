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
