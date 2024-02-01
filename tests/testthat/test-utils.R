test_that("weights_series returns vector of normalised weights", {
  s_params <- list(mean=5, sd=0.01)
  t_max <- 10
  w <- weights_series(t_max, s_params)
  expect_equal(length(w), t_max)
  expect_equal(sum(w), 1)
  expect(which.max(w), 5)
})

test_that("expected_cases returns reasonable values", {
  s_params <- list(mean=5, sd=0.01)
  t_max <- 20
  w <- weights_series(t_max, s_params)
  Rt <- 1
  cases_history <- rep(1, t_max)
  expect_equal(expected_cases(Rt, w, cases_history), 1)

  cases_history <- rep(1, t_max + 2)
  expect_error(expected_cases(Rt, w, cases_history))
})

test_that("pgamma_mean_sd produces densities peaked near correct vals", {
  mu <- 5
  sd <- 0.1
  x <- seq(0, 20, 1)
  max_ind <- which.max(diff(pgamma_mean_sd(x, mu, sd)))
  expect_equal(max_ind, 5)
})

test_that("qgamma_mean_sd produces reasonable inverse", {
  mu <- 5
  sd <- 0.1
  peak <- round(qgamma_mean_sd(0.5, mu, sd))
  expect_equal(peak, 5)
})

test_that("dgamma_mean_sd produces reasonable values", {
  mu <- 5
  sd <- 0.1
  x <- seq(1, 20, 1)
  max_ind <- which.max(dgamma_mean_sd(x, mu, sd))
  expect_equal(max_ind, 5)

  x <- 3
  val <- dgamma_mean_sd(x, mu, sd)
  log_val <- dgamma_mean_sd(x, mu, sd, log=TRUE)
  expect_equal(log(val), log_val)
})

test_that("gamma_discrete_pmf sums near to 1", {
  x <- seq(0, 10000, 1)
  reporting_parameters <- list(mean=1, sd=1)
  vals <- purrr::map_dbl(x, ~gamma_discrete_pmf(., reporting_parameters))
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

test_that("check_parameter_names works as desired", {
  r_params <- list(mean=3, sd=1)
  s_params <- list(mean1=3, sd=1)
  expect_error(check_parameter_names(r_params, s_params))
  r_params <- list(3, 1)
  s_params <- list(mean=3, sd=1)
  expect_error(check_parameter_names(r_params, s_params))
  r_params <- list(mean=3, sd=1)
  s_params <- list(mean=3, sd=1)
  expect_equal(check_parameter_names(r_params, s_params), NULL)
})

