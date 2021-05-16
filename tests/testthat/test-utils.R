test_that("weights_series returns vector of normalised weights", {
  s_params <- list(mean=5, sd=0.01)
  t_max <- 10
  w <- weights_series(t_max, s_params)
  expect_equal(length(w), t_max)
  expect_equal(sum(w), 1)
  expect(which.max(w), 5)
})
