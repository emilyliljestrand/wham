test_that("aic calculates AIC from a WHAM model fit", {
  mod <- list(opt = list(obj = 12.5, par = c(alpha = 1, beta = 2, gamma = 3)))

  expect_equal(aic(mod), 31)
})