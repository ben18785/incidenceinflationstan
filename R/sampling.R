#' Samples a case count arising on a given onset day
#'
#' The distribution being drawn from is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param max_cases maximum possible cases thought to arise on a given day
#' @inheritParams conditional_cases_logp
#' @param ndraws number of draws of cases
#'
#' @return a sampled case count arising on a given onset day
sample_true_cases_single_onset <- function(
  observation_df, cases_history, max_cases,
  Rt, day_onset, serial_parameters, reporting_parameters,
  ndraws=1) {
  max_observed_cases <- max(observation_df$cases_reported)
  if(max_observed_cases > max_cases)
    stop("Max possible cases should be (much) greater than max observed cases.")
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_df, cases_history,
                                  Rt, day_onset, serial_parameters, reporting_parameters)
  probs <- exp(logps - matrixStats::logSumExp(logps))
  sample(possible_cases, ndraws, prob=probs, replace=TRUE)
}

#' Calculates max number of days we are uncertain about reporting
#'
#' @param p_gamma_cutoff a p value (0 <= p <= 1) indicating the threshold above which
#' we deem certainty
#' @inheritParams conditional_cases_logp
#'
#' @return a number of days
max_uncertain_days <- function(p_gamma_cutoff, reporting_parameters) {
  r_mean <- reporting_parameters$mean
  r_sd <- reporting_parameters$sd
  days_from_end <- qgamma(p_gamma_cutoff, r_mean^2 / r_sd^2, r_mean / r_sd^2)
  days_from_end
}

#' Draws a possible history (or histories) of cases
#'
#' The distribution being drawn from at each time t is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param observation_onset_df a tibble with three columns: time_onset, time_reported, cases_reported
#' @inheritParams sample_true_cases_single_onset
#' @inheritParams true_cases
#' @inheritParams max_uncertain_days
#'
#' @return a vector of cases corresponding to the number of cases arising each onset day
#' @export
sample_cases_history <- function(
  observation_onset_df, max_cases,
  Rt_function, serial_parameters, reporting_parameters,
  ndraws=1,
  p_gamma_cutoff=0.99) {

  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_parameters)
  start_uncertain_period <- max(observation_onset_df$time_onset) - uncertain_period
  observation_history_df <- observation_onset_df %>%
    dplyr::group_by(time_onset) %>%
    dplyr::mutate(cases_true=ifelse(time_onset < start_uncertain_period,
                                    max(cases_reported), NA)) %>%
    ungroup()

  onset_times <- unique(observation_history_df$time_onset)
  onset_times_uncertain_period <- onset_times[onset_times >= start_uncertain_period]

  for(i in seq_along(onset_times_uncertain_period)) {
    onset_time <- onset_times_uncertain_period[i]
    observation_df <- observation_history_df %>%
      dplyr::filter(time_onset==onset_time)
    pre_observation_df <- observation_history_df %>%
      dplyr::filter(time_onset < onset_time) %>%
      dplyr::select(time_onset, cases_true) %>%
      unique()
    cases_history <- rev(pre_observation_df$cases_true)
    observation_df <- observation_df %>%
      dplyr::select(time_reported, cases_reported)
    Rt <- Rt_function(onset_time)
    case <- sample_true_cases_single_onset(
      observation_df=observation_df,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=Rt,
      day_onset=onset_time,
      serial_parameters=serial_parameters,
      reporting_parameters=reporting_parameters,
      ndraws=ndraws)
    index_onset_time <- which(observation_history_df$time_onset==onset_time)
    observation_history_df$cases_true[index_onset_time] <- case
  }
  observation_history_df
}
