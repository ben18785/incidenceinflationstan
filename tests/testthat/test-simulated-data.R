test_that("case_undetected_prob returns correct probability", {
  expect_true(abs(case_undetected_prob(20, 10, 10, 5)
                  - 0.4334701) < 0.01)
})
