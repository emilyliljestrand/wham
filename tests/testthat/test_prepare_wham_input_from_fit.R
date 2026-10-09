context("Prepare WHAM input from fitted model")

test_that("input is rebuilt from fitted model options", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  fit_RDS <- tempfile(fileext = ".rds")
  on.exit(unlink(fit_RDS), add = TRUE)
  saveRDS(list(input = input), fit_RDS)

  expect_equal(
    prepare_wham_input_from_fit(fit_RDS)$model_name,
    "WHAM for unnamed stock"
  )
})

test_that("invalid fitted model files produce an informative error", {
  fit_RDS <- tempfile(fileext = ".rds")
  on.exit(unlink(fit_RDS), add = TRUE)
  saveRDS(list(), fit_RDS)

  expect_error(
    prepare_wham_input_from_fit(fit_RDS),
    "does not contain a fitted WHAM model"
  )
})

test_that("models without options return their stored input with a warning", {
  input <- list(data = list(), model_name = "old")
  fit_RDS <- tempfile(fileext = ".rds")
  on.exit(unlink(fit_RDS), add = TRUE)
  saveRDS(list(input = input), fit_RDS)

  expect_warning(
    expect_identical(prepare_wham_input_from_fit(fit_RDS), input),
    "does not contain preparation options"
  )
})
