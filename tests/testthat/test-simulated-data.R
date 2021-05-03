test_that("undetected_prob returns correct probability", {
  expect_true(abs(undetected_prob(20, 10, list(mean=10, sd=5))
                  - 0.4334701) < 0.01)
})

test_that("undetected_prob approaches 0 as t1 approaches inf", {
  expect_equal(undetected_prob(Inf, 10, list(mean=10, sd=5)), 0)
})

test_that("undetected_prob only relative time differences matter", {
  expect_equal(undetected_prob(20, 10, list(mean=10, sd=5)),
               undetected_prob(30, 20, list(mean=10, sd=5)))
})

test_that("undetected_prob increases in mean increase prob", {
  expect_true(undetected_prob(20, 10, list(mean=10, sd=5)) >
                undetected_prob(20, 10, list(mean=5, sd=5)))
})

test_that("undetected_prob throws error if measurement before onset", {
  expect_error(undetected_prob(9, 10, list(mean=10, sd=5)))
})

test_that("detected_after_unobserved_prob returns 0 or 1 as special cases", {
  # d2 = d1
  expect_equal(detected_after_unobserved_prob(2, 2, 0, list(mean=5, sd=5)), 0)
  # d2 = inf
  expect_equal(detected_after_unobserved_prob(Inf, 2, 0, list(mean=5, sd=5)), 1)
})

test_that("detected_after_unobserved_prob increases with latter observation
          time", {
  times <- seq(2, 10, 1)
  probs <- purrr::map_dbl(times,
                          ~detected_after_unobserved_prob(., 2, 1,
                                                          list(mean=5, sd=5)))
  expect_equal(sum(diff(probs) < 0), 0)
})

test_that("detected_after_unobserved_prob throws error if second measurement day
          before first", {
  expect_error(detected_after_unobserved_prob(3, 4, 0, list(mean=5, sd=5)))
})

test_that("observed_cases_single returns a number between 0 and I_remaining", {
  ntries <- 20
  count_fails <- 0
  I_true <- 10
  I_obs <- 2
  I_remaining <- I_true - I_obs
  for(i in 1:ntries) {
    cases_obs <- observed_cases_single(I_obs, I_true, 20, 1, 0,
                                       list(mean=10, sd=5))
    if(cases_obs > I_remaining | cases_obs < 0)
      count_fails <- count_fails + 1
  }
  expect_equal(count_fails, 0)
})

test_that("observed_cases_single throws error if observed cases exceeds true", {
  expect_error(observed_cases_single(15, 13, 23, 22, 20, list(mean=10, sd=5)))
})

test_that("observed_cases_trajectory produces trajectories with given
          characteristics", {
  days_reporting <- seq(1, 10, 1)
  true_cases <- 30
  delay_parameters <- list(mean=1, sd=1)
  day_onset <- 0
  reported_cases <- observed_cases_trajectory(true_cases, days_reporting,
                                              day_onset,
                                              delay_parameters)
  expect_equal(length(reported_cases), length(days_reporting))

  nreps <- 20
  count_increasing <- 0
  count_max <- 0
  for(i in 1:nreps) {
    cases <- observed_cases_trajectory(true_cases, days_reporting, day_onset,
                                       delay_parameters)
    count_increasing <- count_increasing + dplyr::if_else(
      sum(diff(cases) < 0) > 0, 1, 0)
    count_max <- count_max + dplyr::if_else(max(cases) > true_cases, 1, 0)
  }
  expect_equal(count_increasing, 0)
  expect_equal(count_max, 0)
})

test_that("observed_cases_trajectory increases faster with shorter delay", {
  days_reporting <- seq(1, 10, 1)
  true_cases <- 30
  delay_parameters_fast <- list(mean=1, sd=1)
  delay_parameters_slow <- list(mean=40, sd=1)
  day_onset <- 0

  reported_cases_fast <- observed_cases_trajectory(
    true_cases, days_reporting, day_onset, delay_parameters_fast)
  reported_cases_slow <- observed_cases_trajectory(
    true_cases, days_reporting, day_onset, delay_parameters_slow)
  expect_true(max(reported_cases_fast) > max(reported_cases_slow))
})

test_that("gamma_discrete_pmf sums near to 1", {
  x <- seq(0, 10000, 1)
  delay_parameters <- list(mean=1, sd=1)
  vals <- purrr::map_dbl(x, ~gamma_discrete_pmf(., delay_parameters))
  expect_true(1 - sum(vals) < 0.01)
})

test_that("gamma_discrete_pmf peaks where it should", {
  d_params <- list(mean=5, sd=1)
  expect_true(gamma_discrete_pmf(5, d_params) >
              gamma_discrete_pmf(10, d_params))
  d_params <- list(mean=10, sd=1)
  expect_true(gamma_discrete_pmf(5, d_params) <
                gamma_discrete_pmf(10, d_params))
})
