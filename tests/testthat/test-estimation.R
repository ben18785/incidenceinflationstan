test_that("observation_process_logp returns reasonable log-probs", {
  reporting_parameters <- list(mean=5, sd=4)
  day_2 <- 4
  day_1 <- 2
  day_onset <- 0
  I_true <- 10
  I_day_1 <- 3

  # I_day_2 > I_true should be impossible
  I_day_2 <- 11
  logp <- observation_process_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 < I_day_1 should be impossible
  I_day_2 <- 2
  logp <- observation_process_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 near I_true should be more likely with shorter delays
  I_day_2 <- 10
  logp <- observation_process_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  reporting_parameters1 <- list(mean=1, sd=1)
  logp1 <- observation_process_logp(I_true, I_day_2, I_day_1,
                                    day_2, day_1, day_onset, reporting_parameters1)
  expect_true(logp1 > logp)

  # prob 1 if I_day_1=I_true
  I_day_2 <- 10
  I_day_1 <- 10
  logp <- observation_process_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(0, logp)
})

test_that("state_process_logp produces reasonable log prob values", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), -1)

  Rt <- 2
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 2
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), log(2 * exp(-2)))
})
