ecov_safety_input <- function(years = 2000:2005) {
  list(
    data = list(
      n_stocks = 1L,
      n_indices = 1L,
      n_ages = 2L,
      n_regions = 1L,
      n_seasons = 1L,
      n_fleets = 1L,
      n_years_model = length(years),
      recruit_model = 2L,
      NAA_re_model = 1L,
      FXSPR_init = rep(0.1, length(years)),
      FMSY_init = rep(0.1, length(years))
    ),
    par = list(),
    map = list(),
    log = list(),
    options = list(),
    years = years,
    by_pwi = TRUE
  )
}

ecov_safety_options <- function(
  years = 2000:2005,
  n_ecov = 1L,
  logsigma = "est_re"
) {
  list(
    label = paste0("covariate ", seq_len(n_ecov)),
    mean = matrix(seq_len(length(years) * n_ecov), length(years), n_ecov),
    use_obs = matrix(1, length(years), n_ecov),
    year = years,
    logsigma = logsigma,
    process_model = "ar1"
  )
}

test_that("est_re sigma parameters have finite starts and independent annual maps", {
  for (n_ecov in c(1L, 3L)) {
    input <- set_ecov(
      ecov_safety_input(2001:2005),
      ecov_safety_options(n_ecov = n_ecov)
    )
    expect_equal(dim(input$par$Ecov_obs_sigma_par), c(2L, n_ecov))
    expect_equal(nlevels(input$map$Ecov_obs_sigma_par), 2L * n_ecov)
    expect_equal(dim(input$par$Ecov_obs_logsigma_re), c(6L, n_ecov))
    expect_equal(nlevels(input$map$Ecov_obs_logsigma_re), 6L * n_ecov)
    expect_equal(input$par$Ecov_obs_logsigma_re, matrix(-1.3, 6L, n_ecov))
    expect_equal(all(is.finite(input$par$Ecov_obs_logsigma)), TRUE)
    expect_equal(all(is.finite(input$par$Ecov_obs_sigma_par)), TRUE)
    expect_equal("Ecov_obs_logsigma_re" %in% input$random, TRUE)
  }
})

test_that("sigma RE dimensions and maps include both ends of padded years", {
  ecov <- ecov_safety_options(years = 2002:2003, n_ecov = 2L)
  ecov$logsigma <- list(
    matrix(c(-0.7, -0.9, NA, -1.1), 2L, 2L),
    c("est_re", "est_re")
  )
  input <- set_ecov(ecov_safety_input(), ecov)
  expect_equal(input$years_Ecov, 1999:2005)
  expect_equal(dim(input$par$Ecov_obs_logsigma_re), c(7L, 2L))
  expect_length(input$map$Ecov_obs_logsigma_re, 14L)
  expect_equal(nlevels(input$map$Ecov_obs_logsigma_re), 14L)
  expect_equal(
    input$par$Ecov_obs_logsigma_re[, 1],
    c(-0.7, -0.7, -0.7, -0.7, -0.9, -0.9, -0.9)
  )
  expect_equal(all(is.finite(input$par$Ecov_obs_logsigma_re)), TRUE)
  expect_equal(colSums(input$data$Ecov_use_obs), c(2, 2))
})

test_that("mixed fixed and estimated sigma options survive a rebuild", {
  fixed <- matrix(seq(-2.4, -0.1, length.out = 24L), 6L, 4L)
  ecov <- ecov_safety_options(
    n_ecov = 4L,
    logsigma = list(fixed, c(NA, "est_1", "est_re", NA))
  )
  first <- set_ecov(ecov_safety_input(), ecov)
  second <- set_ecov(first, first$options$ecov)
  expect_identical(first$options$ecov$logsigma, ecov$logsigma)
  expect_identical(second$options$ecov$logsigma, ecov$logsigma)
  expect_equal(first$par$Ecov_obs_logsigma[-1, c(1, 4)], fixed[, c(1, 4)])
  expect_identical(second$par$Ecov_obs_logsigma, first$par$Ecov_obs_logsigma)
  expect_identical(
    second$par$Ecov_obs_logsigma_re,
    first$par$Ecov_obs_logsigma_re
  )
  expect_identical(
    second$map$Ecov_obs_logsigma_re,
    first$map$Ecov_obs_logsigma_re
  )
  expect_equal(first$data$Ecov_obs_sigma_opt, c(1, 2, 4, 1))
  expect_equal(nlevels(first$map$Ecov_obs_sigma_par), 2L)
  expect_equal(nlevels(first$map$Ecov_obs_logsigma_re), 7L)
})

test_that("fixed sigma configurations allocate all covariate columns", {
  for (logsigma in list(c(-1, -2), matrix(c(-1, -2), 1L, 2L))) {
    input <- set_ecov(
      ecov_safety_input(),
      ecov_safety_options(n_ecov = 2L, logsigma = logsigma)
    )
    expect_equal(dim(input$par$Ecov_obs_sigma_par), c(2L, 2L))
    expect_equal(dim(input$par$Ecov_obs_logsigma_re), c(7L, 2L))
    expect_equal(
      input$par$Ecov_obs_logsigma,
      matrix(rep(c(-1, -2), each = 7L), 7L, 2L)
    )
    expect_equal(nlevels(input$map$Ecov_obs_sigma_par), 0L)
    expect_equal(nlevels(input$map$Ecov_obs_logsigma_re), 0L)
  }
})

test_that("invalid environmental years fail before padding", {
  invalid <- list(
    numeric(),
    c(2000, NA),
    c(2000, Inf),
    c(2000, NaN),
    c(2000, 2000),
    c(2001, 2000),
    c(2000, 2002),
    c("2000", "2001"),
    matrix(2000:2001, 2L)
  )
  for (years in invalid) {
    ecov <- ecov_safety_options()
    ecov$year <- years
    error <- tryCatch(set_ecov(ecov_safety_input(), ecov), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "nonempty, finite numeric vector of consecutive years",
      fixed = TRUE
    )
  }
  ecov <- ecov_safety_options()
  ecov$year <- 2000:2001
  error <- tryCatch(set_ecov(ecov_safety_input(), ecov), error = identity)
  expect_match(conditionMessage(error), "not the same length", fixed = TRUE)
})
