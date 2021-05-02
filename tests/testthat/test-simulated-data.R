test_that("undetected_prob returns correct probability", {
  expect_true(abs(undetected_prob(20, 10, 10, 5)
                  - 0.4334701) < 0.01)
})

test_that("undetected_prob approaches 0 as t1 approaches inf", {
  expect_equal(undetected_prob(Inf, 10, 10, 5), 0)
})

test_that("undetected_prob only relative time differences matter", {
  expect_equal(undetected_prob(20, 10, 10, 5),
               undetected_prob(30, 20, 10, 5))
})

test_that("undetected_prob increases in mean increase prob", {
  expect_true(undetected_prob(20, 10, 10, 5) > undetected_prob(20, 10, 5, 5))
})

test_that("detected_after_unobserved_prob returns 0 or 1 as special cases", {
  # d2 = d1
  expect_equal(detected_after_unobserved_prob(2, 2, 0, 10, 5), 0)
  # d2 = inf
  expect_equal(detected_after_unobserved_prob(Inf, 2, 0, 10, 5), 1)
})

test_that("detected_after_unobserved_prob increases with latter observation
          time", {
  times <- seq(2, 10, 1)
  probs <- purrr::map_dbl(times,
                          ~detected_after_unobserved_prob(., 2, 1, 10, 5))
  expect_equal(sum(diff(probs) < 0), 0)
})
