test_that("projection sigma extension leaves inactive inputs unchanged", {
  par <- list(Ecov_obs_logsigma_re = matrix(3, 2, 1))
  map <- list(Ecov_obs_logsigma_re = factor(c(NA, NA)))
  data <- list(Ecov_obs_sigma_opt = 1L, n_Ecov = 1L)

  expect_identical(
    extend_ecov_input(par, map, data, 2L),
    list(par = par, map = map)
  )
})

test_that("projection sigma extension preserves existing sharing", {
  par <- list(
    Ecov_obs_logsigma_re = matrix(c(1, 2, 3, 4), 2, 2),
    Ecov_obs_sigma_par = matrix(c(-1, -2, -3, -4), 2, 2)
  )
  map <- list(Ecov_obs_logsigma_re = factor(c(1, 1, NA, NA)))
  data <- list(Ecov_obs_sigma_opt = c(4L, 1L), n_Ecov = 2L)

  extended <- extend_ecov_input(par, map, data, 2L)
  mapping <- matrix(as.character(extended$map$Ecov_obs_logsigma_re), 4, 2)

  expect_identical(
    extended$par$Ecov_obs_logsigma_re[1:2, ],
    par$Ecov_obs_logsigma_re
  )
  expect_equal(
    extended$par$Ecov_obs_logsigma_re[3:4, ],
    matrix(c(-1, -1, 0, 0), 2, 2)
  )
  expect_identical(mapping[1:2, 1], c("1", "1"))
  expect_identical(mapping[3:4, 1], c("2", "3"))
  expect_identical(mapping[, 2], rep(NA_character_, 4))
})

test_that("report name lookup prefers the current name", {
  report <- matrix(
    0,
    3,
    1,
    dimnames = list(c("log_SSB", "log_SSB_all", "other"), NULL)
  )

  expect_identical(wham_report_index(report, "log_SSB_all", "log_SSB"), 2L)
  expect_identical(wham_report_index(report, "absent", "log_SSB"), 1L)
  expect_identical(wham_report_index(report, "absent", "missing"), integer())
})
