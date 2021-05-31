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
  if(dplyr::last(probs) > 0.01)
    warning(paste0("Cases too few for onset day: ", day_onset,
                   ". Increase max_cases."))
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
  days_from_end <- qgamma_mean_sd(p_gamma_cutoff, r_mean, r_sd)
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
#' @return a tibble with an extra cases_estimated column
#' @export
#' @importFrom rlang .data
sample_cases_history <- function(
  observation_onset_df, max_cases,
  Rt_function, serial_parameters, reporting_parameters,
  p_gamma_cutoff=0.99) {

  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_parameters)
  start_uncertain_period <- max(observation_onset_df$time_onset) - uncertain_period
  observation_history_df <- observation_onset_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::mutate(cases_estimated=ifelse(.data$time_onset < start_uncertain_period,
                                    max(.data$cases_reported), NA)) %>%
    dplyr::ungroup()
  onset_times <- unique(observation_history_df$time_onset)
  onset_times_uncertain_period <- onset_times[onset_times >= start_uncertain_period]

  for(i in seq_along(onset_times_uncertain_period)) {
    onset_time <- onset_times_uncertain_period[i]
    snapshots_at_onset_time_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset==onset_time) %>%
      dplyr::select(.data$time_reported, .data$cases_reported)
    pre_observation_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::select(.data$time_onset, .data$cases_estimated) %>%
      unique() %>%
      dplyr::arrange(dplyr::desc(.data$time_onset))
    cases_history <- pre_observation_df$cases_estimated
    Rt <- Rt_function(onset_time)
    case <- sample_true_cases_single_onset(
      observation_df=snapshots_at_onset_time_df,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=Rt,
      day_onset=onset_time,
      serial_parameters=serial_parameters,
      reporting_parameters=reporting_parameters,
      ndraws=1)
    index_onset_time <- which(observation_history_df$time_onset==onset_time)
    observation_history_df$cases_estimated[index_onset_time] <- case
  }
  observation_history_df
}

#' Title
#'
#' @param R0_piece_index the index of the R0 piece being sampled
#' @param cases_history_df a tibble with three columns: time_onset, cases_true
#' and R0_index
#' @param R0_prior_parameters a list with elements 'shape' and 'rate' describing
#' the gamma prior for R0
#' @inheritParams sample_cases_history
#' @param ndraws number of draws of R0
#'
#' @return a draw for R0
#' @importFrom rlang .data
sample_R0_single_piece <- function(R0_piece_index,
                                   cases_history_df,
                                   R0_prior_parameters,
                                   serial_parameters,
                                   serial_max=40,
                                   ndraws=1) {
  short_df <- cases_history_df %>%
    dplyr::filter(.data$R0_index <= R0_piece_index)
  time_max_post_initial_period <- max(short_df$time_onset) - serial_max

  # sample from prior since no data
  alpha <- R0_prior_parameters$shape
  beta <- R0_prior_parameters$rate
  if(time_max_post_initial_period <= 0)
    return(stats::rgamma(ndraws, alpha, beta))

  # if some data but not enough for whole period
  # do not use truncated points as observed data
  # (but they will be used as covariates)
  short_df <- short_df %>%
    dplyr::mutate(time_after_start = .data$time_onset - serial_max) %>%
    dplyr::mutate(is_observed_data=if_else(
      (time_after_start > 0) & (R0_index == R0_piece_index), 1, 0))
  onset_times <- short_df %>%
    dplyr::filter(is_observed_data == 1) %>%
    dplyr::pull(time_onset)

  posterior_shape <- alpha
  posterior_rate <- beta
  w <- weights_series(serial_max, serial_parameters)
  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    true_cases <- short_df %>%
      dplyr::filter(.data$time_onset == onset_time) %>%
      dplyr::pull(cases_true)
    posterior_shape <- posterior_shape + true_cases
    cases_history <- short_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::arrange(desc(time_onset)) %>%
      dplyr::pull(cases_true)
    cases_history <- cases_history[1:serial_max]
    posterior_rate <- posterior_rate + sum(w * cases_history)
  }
  stats::rgamma(ndraws, posterior_shape, posterior_rate)
}
