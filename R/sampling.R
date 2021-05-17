#' Samples a case count arising on a given onset day
#'
#' The distribution being drawn from is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param cases_true a vector of possible case counts
#' @inheritParams conditional_cases_logp
#' @param ndraws number of draws of cases
#'
#' @return a sampled case count arising on a given onset day
sample_true_cases <- function(cases_true, observation_df, cases_history,
                              Rt, day_onset, serial_parameters, reporting_parameters,
                              ndraws=1) {
  logps <- conditional_cases_logp(cases_true, observation_df, cases_history,
                                  Rt, day_onset, serial_parameters, reporting_parameters)
  probs <- exp(logps - matrixStats::logSumExp(logps))
  sample(cases_true, ndraws, prob=probs, replace=TRUE)
}
