#' Calculates observation process log probability density for a pair of days
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
#' @param cases_day_2 number of cases arising on day_onset observed by day_2
#' @param cases_day_1 number of cases arising on day_onset observed by day_1 < day_2
#' @param cases_true true case count(s) originating on day_onset: can be a single value or vector
#' @inheritParams detected_after_unobserved_prob
#' @inheritParams observed_cases_single
#'
#' @return a log-probability or vector of log-probabilities
observation_process_single_logp <- function(cases_true, cases_day_2, cases_day_1, day_2, day_1,
                                            day_onset, reporting_parameters){
  cases_observed <- cases_day_2 - cases_day_1
  cases_remaining <- cases_true - cases_day_1
  p_detect <- detected_after_unobserved_prob(day_2, day_1, day_onset,
                                             reporting_parameters)
  stats::dbinom(cases_observed, cases_remaining, p_detect, log = T)
}

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
#' This function calculates the combined probability of all observed cases pertaining
#' to a particular onset date.
#'
#' @param observation_df a two-column data frame with columns 'time_reported' and 'cases_reported'
#' @inheritParams observation_process_single_logp
#'
#' @return a log-probability or vector of log-probabilities
observation_process_logp <- function(observation_df, cases_true,
                                     day_onset, reporting_parameters){
  observation_days <- observation_df$time_reported
  observation_cases <- observation_df$cases_reported
  change_days <- diff(observation_days)
  if(any(change_days <= 0))
    stop("days in observation matrix must be increasing.")
  change_cases <- diff(observation_cases)
  if(any(change_cases < 0))
    stop("reported cases in observation matrix must be increasing.")

  num_periods <- nrow(observation_df) - 1
  logp <- 0
  for(i in 1:num_periods) {
    logp_day <- observation_process_single_logp(cases_true=cases_true,
                                         cases_day_2=observation_cases[i + 1],
                                         cases_day_1=observation_cases[i],
                                         day_2=observation_days[i + 1],
                                         day_1=observation_days[i],
                                         day_onset=day_onset,
                                         reporting_parameters=reporting_parameters)
    logp <- logp + logp_day
  }
  logp
}

#' Determines the probability of a given number of cases given a history of them
#'
#' Calculates:
#' \deqn{Pois(cases_true_t|Rt * \sum_tau=1^t_max w_t cases_true_t-tau)}
#'
#' @inheritParams observed_cases_single
#' @inheritParams true_cases_single
#' @inheritParams gamma_discrete_pmf
#'
#' @return a log-probability or vector of log-probabilities
state_process_logp <- function(cases_true, cases_history, Rt, serial_parameters){
  t_max <- length(cases_history)
  w <- weights_series(t_max, serial_parameters)
  mean_cases <- expected_cases(Rt=Rt, weights=w, cases_history=cases_history)
  stats::dpois(cases_true, lambda=mean_cases, log = T)
}

#' Calculates (unnormalised) log-probability of a true case count pertaining to cases arising
#' on a particular onset day
#'
#' There are two components to this:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @inheritParams state_process_logp
#' @inheritParams observation_process_logp
#'
#' @return an (unnormalised) log-probability or vector of such log-probabilities
conditional_cases_logp <- function(cases_true, observation_df, cases_history,
                                   Rt, day_onset, serial_parameters, reporting_parameters) {

  logp_observation <- 0
  if(nrow(observation_df) > 1) # handling cases arising today which were reported today
    logp_observation <- observation_process_logp(observation_df=observation_df,
                                                 cases_true=cases_true,
                                                 day_onset=day_onset,
                                                 reporting_parameters=reporting_parameters)
  logp_state <- state_process_logp(cases_true=cases_true,
                                   cases_history=cases_history,
                                   Rt=Rt,
                                   serial_parameters=serial_parameters)
  logp_observation + logp_state
}

#' Observation probability across all onset times
#'
#' @param snapshot_with_true_cases_df a tibble with four columns:
#' time_onset, time_reported, cases_reported, cases_true
#' @param reporting_parameters a tibble with columns: 'reporting_piece_index', 'mean', 'sd'
#'
#' @return a log probability
#' @importFrom rlang .data
observation_process_all_times_logp <- function(
  snapshot_with_true_cases_df,
  reporting_parameters){

  if(methods::is(reporting_parameters, "list"))
    stop("Reporting parameters should be a tibble not a list.")

  if(!"reporting_piece_index" %in% colnames(reporting_parameters))
    stop("reporting_parameters must contain a column: 'reporting_piece_index'.")

  reporting_parameters <- snapshot_with_true_cases_df %>%
    dplyr::left_join(reporting_parameters, by="reporting_piece_index") %>%
    dplyr::select("time_onset", "mean", "sd") %>%
    unique()

  logp <- 0
  onset_times <- unique(snapshot_with_true_cases_df$time_onset)

  # check that only one true case per onset time
  n_onset <- length(onset_times)
  test_df <- snapshot_with_true_cases_df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique()
  if(n_onset != nrow(test_df))
    stop("There must be only one true case measurement per each onset time.")

  mu <- reporting_parameters$mean
  sd <- reporting_parameters$sd
  if(sum(mu < 0) > 0 | sum(sd < 0) > 0)
    return(-Inf)

  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    short_df <- snapshot_with_true_cases_df %>%
      dplyr::filter(.data$time_onset==onset_time)
    cases_true <- short_df$cases_true[1]
    reporting_parameters_current <- list(mean=reporting_parameters$mean[i],
                                         sd=reporting_parameters$sd[i])
    if(nrow(short_df) > 1) {# handling cases arising today which were reported today
      logp_new <- observation_process_logp(
        short_df,
        cases_true=cases_true,
        day_onset=onset_time,
        reporting_parameters=reporting_parameters_current)
      logp <- logp + logp_new
    }
  }
  logp
}
