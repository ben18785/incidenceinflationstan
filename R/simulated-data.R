#' Probability case remains undetected over time
#'
#' @param t_primed time subsequent to time when case originates
#' @param t time when case originates
#' @param delay_mean mean of gamma distribution governing reporting delay distribution
#' @param delay_sd sd of gamma distribution governing reporting delay distribution
#'
#' @return a probability
case_undetected_prob <- function(t_primed, t,
                                 delay_mean, delay_sd){
  1 - stats::pgamma(t_primed - t, delay_mean^2 / delay_sd^2,
                    delay_mean / delay_sd^2)
}
