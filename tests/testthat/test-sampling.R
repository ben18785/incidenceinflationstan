test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 10
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  w <- weights_series(t_max, s_params)
  mean_cases <- expected_cases(Rt, w, cases_history)
  case <- sample_true_cases_single_onset(observation_df=observation_matrix,
                                 cases_history=cases_history,
                                 max_cases=max_cases,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(case >= min(observation_matrix$cases_reported))
  expect_true(case <= max_cases)

  # when reporting mean is low, expect to have pretty much seen all cases
  reporting_parameters <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                            cases_history=cases_history,
                            max_cases=max_cases,
                            Rt=Rt,
                            day_onset=day_onset,
                            serial_parameters=s_params,
                            reporting_parameters=reporting_parameters,
                            ndraws=20)
  expect_true(any(cases <= 5))

  # check that cases jump if serial interval looks further back
  cases_history <- c(rep(4, t_max / 2), rep(1000, t_max / 2))
  s_params <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  s_params <- list(mean=25, sd=1)
  cases1 <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  expect_true(mean(cases1) > mean(cases))
})
