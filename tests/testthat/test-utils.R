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
