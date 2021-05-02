#' Probability case remains undetected over time
#'
#' This is effectively a survival function, which shows the
#' probability that an actual case originating on day t0
#' remains unreported at t1.
#'
#' @param t1 time subsequent to time when case originates, t1>=t0
#' @param t0 time when case originates
#' @param delay_mean mean of gamma distribution governing reporting delay distribution
#' @param delay_sd std. dev. of gamma distribution governing reporting delay distribution
#'
#' @return a probability
undetected_prob <- function(t1, t0, delay_mean, delay_sd){
  1 - stats::pgamma(t1 - t0, delay_mean^2 / delay_sd^2,
                    delay_mean / delay_sd^2)
}

#' Probability a case is detected by a later time given it was unobserved
#' earlier
#'
#' This is the probability that a case originating on day d0 is detected between
#' days d1 and d2 given that it was unobserved on day d1. This probability is
#' calculated using the survival functions (which indicate the probability a
#' case remains undetected by a certain time):
#' \deqn{prob(detect) = \frac{S(d1|d0) - S(d2|d0)}{S(d1|d0)}}
#'
#' @param d2 day of last observation d2 > d1 >= d0
#' @param d1 day of first observation d1 >= d0
#' @param d0 day when case originates
#' @inheritParams undetected_prob
#'
#' @return a probability
detected_after_unobserved_prob <- function(d2, d1, d0, delay_mean, delay_sd){
  (undetected_prob(d1 + 0.5, d0, delay_mean, delay_sd) -
   undetected_prob(d2 + 0.5, d0, delay_mean, delay_sd)) /
  undetected_prob(d1 + 0.5, d0, delay_mean, delay_sd)
}
