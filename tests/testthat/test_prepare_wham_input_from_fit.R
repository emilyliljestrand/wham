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

test_that("missing recruitment decoupling option defaults to coupled", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  input$options$NAA_re$decouple_recruitment <- NULL

  rebuilt <- prepare_wham_input_from_fit(list(input = input))

  expect_identical(rebuilt$data$decouple_recruitment, 0L)
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

test_that("refitting supplies missing year indices and preserves stored indices", {
  fit <- list(
    input = list(data = list(n_years_model = 3L)),
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
  )
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) {
      result <- fit
      result$input <- input
      result
    }
  ), parent = environment(refit_wham))

  expect_warning(
    capture.output(result <- refit(fit)),
    "does not contain preparation options"
  )
  expect_identical(result$input$data$years_use, 0:2)
  expect_identical(result$input$data$decouple_recruitment, 0L)
  expect_null(fit$input$data$years_use)

  fit$input$data$years_use <- c(0L, 2L)
  expect_warning(
    capture.output(result <- refit(fit)),
    "does not contain preparation options"
  )
  expect_identical(result$input, fit$input)
})

test_that("unsupported legacy configurations are rejected before rebuilding", {
  input <- list(
    asap3 = list(n_ages = 5L),
    data = list(
      Ecov_model = 0, selblock_models_re = 1, M_re_model = 1,
      M_model = 2, decouple_rec = 0, use_steepness = 0, do_proj = 0
    )
  )
  unsupported <- list(
    Ecov_model = 1, selblock_models_re = 2, M_re_model = 2,
    M_model = 3, decouple_rec = 1, use_steepness = 1, do_proj = 1
  )
  for (field in names(unsupported)) {
    legacy_input <- input
    legacy_input$data[[field]] <- unsupported[[field]]
    expect_error(
      migrate_legacy_refit_input(list(input = legacy_input)),
      "Cannot safely migrate this legacy WHAM input"
    )
  }
})

test_that("legacy Butterfish migration preserves the likelihood and estimates", {
  fit_file <- test_path("..", "..", "oldfitobjects", "Butterfishfit.rds")
  skip_if_not(file.exists(fit_file), "Local legacy fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  input <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(fit$input, original_input)
  expect_equal(input$data$years_use, 0:34)
  expect_equal(input$data$bias_correct_pe, fit$input$data$bias_correct_pe)
  expect_equal(input$data$bias_correct_oe, fit$input$data$bias_correct_oe)
  expect_identical(input$data$obsvec, fit$input$data$obsvec)
  expect_identical(input$map$logit_q, fit$input$map$logit_q)
  expect_equal(as.vector(input$par$Mpars), fit$parList$M_a)
  expect_equal(as.vector(input$par$trans_NAA_rho), 2 * fit$parList$trans_NAA_rho)

  model <- TMB::MakeADFun(
    input$data, input$par, DLL = "wham", random = input$random,
    map = input$map, silent = TRUE
  )
  expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)
  report <- model$report()
  for (field in c("nll_NAA", "nll_agg_catch", "nll_agg_indices", "nll_catch_acomp",
                  "nll_index_acomp", "SSB", "NAA", "pred_catch", "pred_indices")) {
    expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]), tolerance = 1e-8)
  }

  fit$call <- quote(fit_wham(input = input, do.retro = FALSE, do.osa = FALSE,
    do.sdrep = FALSE, do.brps = FALSE, do.post.samp = FALSE, MakeADFun.silent = TRUE))
  output <- capture.output(result <- suppressMessages(refit_wham(fit)))
  expect_equal(result$opt$convergence, 0)
  expect_equal(as.numeric(result$opt$objective), as.numeric(fit$opt$objective), tolerance = 1e-8)
  expect_true(any(grepl("max_abs_par_diff", output)))
})
