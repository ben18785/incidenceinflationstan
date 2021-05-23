test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
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

  max_cases <- 2
  expect_warning(sample_true_cases_single_onset(observation_df=observation_matrix,
                                                cases_history=cases_history,
                                                max_cases=max_cases,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters,
                                                ndraws=20))
  max_cases <- 0
  expect_error(sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    ndraws=20))
})

test_that("max_uncertain_days selects reasonable gamma distribution quanties", {
  r_params <- list(mean=7, sd=0.01)
  expect_equal(round(max_uncertain_days(0.5, r_params)), 7)
  expect_true(max_uncertain_days(0.95, r_params) > 7)

  r_params <- list(mean=7, sd=5)
  expect_true(max_uncertain_days(0.95, r_params) > 10)
  expect_equal(max_uncertain_days(0, r_params), 0)
})

test_that("sample_cases_history adds cases_estimated that look reasonable", {
  days_total <- 100
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T)
  max_cases <- 5000
  df_est <- sample_cases_history(df, max_cases, Rt_function, s_params, r_params)
  expect_true(dplyr::all_equal(df_est %>% dplyr::select(-cases_estimated),
                               df))
  expect_true(max(df_est$cases_estimated) < max_cases)
  expect_true(min(df_est$cases_estimated) > 0)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  df_group <- df_est %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(
      cases_reported=max(cases_reported),
      cases_estimated=min(cases_estimated),
      cases_true=mean(cases_true)) %>%
    dplyr::mutate(diff=cases_estimated-cases_reported) %>%
    dplyr::mutate(diff_true=cases_estimated-cases_true)
  expect_equal(sum(df_group$diff < 0), 0)
  expect_true(max(abs(df_group$diff_true)) < 300)
})
