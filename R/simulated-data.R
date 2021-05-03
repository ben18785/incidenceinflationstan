#' Probability case remains undetected over time
#'
#' This is effectively a survival function, which shows the
#' probability that an actual case originating on day t_onset
#' remains unreported at t_1.
#'
#' @param t_1 observation time subsequent to time when case originates: t1>=t0
#' @param t_onset time when case originates
#' @param delay_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the reporting delay distribution
#'
#' @return a probability
#' @examples
#' incidenceinflation:::undetected_prob(2, 1, list(mean=10, sd=5))
undetected_prob <- function(t_1, t_onset, delay_parameters){
  if(t_1 < t_onset)
    stop("t_1 must equal or exceed t_onset")
  delay_mean <- delay_parameters$mean
  delay_sd <- delay_parameters$sd
  1 - stats::pgamma(t_1 - t_onset, delay_mean^2 / delay_sd^2,
                    delay_mean / delay_sd^2)
}

#' Probability a case is detected by a later time given it was unobserved
#' earlier
#'
#' This is the probability that a case originating on day day_onset is detected
#' between days day_1 and day_2 given that it was unobserved on day day_1. This
#' probability is calculated using the survival functions (which indicate the
#' probability a case remains undetected by a certain time):
#' \deqn{prob(detect) = \frac{S(day_1|day_onset) - S(day_2|day_onset)}{S(day_1|day_onset)}}
#'
#' @param day_2 day of last observation: day_2 > day_1 >= day_onset
#' @param day_1 day of first observation: day_1 >= day_onset
#' @param day_onset day when case originates
#' @inheritParams undetected_prob
#'
#' @return a probability
detected_after_unobserved_prob <- function(day_2, day_1, day_onset,
                                           delay_parameters){
  if(day_2 < day_1)
    stop("second observation date must be at same time or after first")
  (undetected_prob(day_1 + 0.5, day_onset, delay_parameters) -
   undetected_prob(day_2 + 0.5, day_onset, delay_parameters)) /
  undetected_prob(day_1 + 0.5, day_onset, delay_parameters)
}

#' Cases arising on a given day which are reported between two subsequent days
#'
#' The number of cases originating on day_onset which are subsequently reported
#' between day_1 and day_2 is modelled as a binomial distribution:
#' \deqn{I_reported\sim binomial(I_true-I_observed, p_detect)}
#' where I_observed is the number of cases originating on day_onset reported on
#' day_1; p_detect is the probability a case was undetected on day_1 but
#' is detected by day_2.
#'
#' @inheritParams detected_after_unobserved_prob
#' @param I_true true case count originating on day_onset
#' @param I_observed observed case count for cases originating on day_onset on
#' day_1
#'
#' @return a count representing number of cases reported between day_1 and day_2
observed_cases_single <- function(I_observed, I_true, day_2, day_1, day_onset,
                                  delay_parameters){
  if(I_true < I_observed)
    stop("true case count must exceed reported")
  I_remaining <- I_true - I_observed
  stats::rbinom(1, I_remaining,
                detected_after_unobserved_prob(day_2, day_1, day_onset,
                                               delay_parameters))
}


#' Generate trajectory of reported cases for a given count originating on a
#' particular day
#'
#' @inheritParams observed_cases_single
#' @param days_reporting a vector of days on which cases were counted
#'
#' @return a vector of reported case counts of same length as days_reporting
observed_cases_trajectory <- function(I_true, days_reporting, day_onset,
                                      delay_parameters){

  I_observed <- vector(length = length(days_reporting))
  I_observed[1] <- 0
  for(i in 1:length(days_reporting)){
    if(i == 1) {
      I_previous_obs <- 0
      day_previous_report <- 0
    } else {
      I_previous_obs <- I_observed[i - 1]
      day_previous_report <- days_reporting[i - 1]
    }
    new_cases <- observed_cases_single(I_previous_obs, I_true,
                                       days_reporting[i],
                                       day_previous_report,
                                       day_onset,
                                       delay_parameters)
    I_observed[i] <- I_previous_obs + new_cases
  }
  I_observed
}

#' Discrete gamma probability mass function
#'
#' @param day day to evaluate pmf
#' @inheritParams undetected_prob
#'
#' @return a probability
gamma_discrete_pmf <- function(day, delay_parameters){
  delay_mean <- delay_parameters$mean
  delay_sd <- delay_parameters$sd
  a <- delay_mean^2 / delay_sd^2
  b <- delay_mean / delay_sd^2
  stats::pgamma(day + 0.5, a, b) - stats::pgamma(day - 0.5, a, b)
}
