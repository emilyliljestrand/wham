context("Prepare WHAM input from fitted model")

test_that("input is rebuilt from fitted model options", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  fit <- list(input = input)

  expect_equal(
    prepare_wham_input_from_fit(fit)$model_name,
    "WHAM for unnamed stock"
  )
})

test_that("invalid fitted model files produce an informative error", {
  expect_error(
    prepare_wham_input_from_fit(list()),
    "does not contain an \\$input component"
  )
})

test_that("models without options return their stored input with a warning", {
  input <- list(data = list(), model_name = "old")

  expect_warning(
    expect_identical(prepare_wham_input_from_fit(list(input = input)), input),
    "does not contain preparation options"
  )
})
