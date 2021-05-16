#' Calculates observation process log probability density
#'
#' The probability of observations is a function of the number of cases remaining
#' to be reported:
#' \deqn{I_remaining = I_true - I_day_1}
#' and the number of cases observed between day_1 and day_2:
#' \deqn{I_obs = I_day_2 - I_day_1}
#' It is given by:
#' \deqn{binomial(I_obs|I_remaining, p_detect)}
#' where p_detect is the probability a thus undetected case is detected between
#' day_1 and day_2.
#'
#'
#' @param cases_day_2 number of cases arising on day_onset observed by day_2
#' @param cases_day_1 number of cases arising on day_onset observed by day_1 < day_2
#' @inheritParams detected_after_unobserved_prob
#' @inheritParams observed_cases_single
#'
#' @return a log-probability
observation_process_logp <- function(cases_true, cases_day_2, cases_day_1, day_2, day_1,
                                           day_onset, reporting_parameters){
  stats::dbinom(cases_day_2 - cases_day_1, cases_true - cases_day_1,
                detected_after_unobserved_prob(day_2, day_1, day_onset,
                                               reporting_parameters),
                log = T)
}


#' Determines the probability of a given number of cases given a history of them
#'
#'
#'
#' @inheritParams observed_cases_single
#' @inheritParams true_cases_single
#' @inheritParams gamma_discrete_pmf
#'
#' @return a log-probability
state_process_logp <- function(cases_true, cases_history, Rt, serial_parameters){
  t_max <- length(cases_history)
  w <- weights_series(t_max, serial_parameters)
  a_sum <- Rt * sum(w * cases_history)
  stats::dpois(cases_true, lambda=a_sum, log = T)
}
