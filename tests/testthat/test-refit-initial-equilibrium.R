test_that("Run4 refits preserve historical equilibrium likelihood and estimates", {
  path <- test_path("..", "..", "oldfitobjects", "BSB_2023Run4fit.RDS")
  skip_if_not(file.exists(path), "Local BSB Run4 fit is unavailable")
  fit <- readRDS(path)
  original_input <- fit$input
  fit$call <- quote(fit_wham(
    input = input,
    do.retro = FALSE,
    do.osa = FALSE,
    do.sdrep = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  ))
  capture.output(result <- suppressMessages(refit_wham(fit)))
  expect_equal(result$input$data$N1_model, -1L)
  expect_equal(
    result$opt$objective,
    fit$opt$objective,
    tolerance = 1e-6,
    scale = 1
  )
  expect_equal(result$opt$par, fit$opt$par, tolerance = 1e-6, scale = 1)
  expect_equal(result$rep$N1, fit$rep$N1, tolerance = 1e-8)
  expect_equal(result$rep$NAA, fit$rep$NAA, tolerance = 1e-8)
  expect_identical(fit$input, original_input)

  input <- result$input
  expect_identical(refit_initial_equilibrium(input, result$rep$N1), input)
  input$data$N1_model <- 1L
  model <- TMB::MakeADFun(
    input$data,
    input$par,
    map = input$map,
    random = input$random,
    DLL = "wham",
    silent = TRUE
  )
  on.exit(TMB::FreeADFun(model))
  model$fn(model$par)
  current_N1 <- model$report()$N1
  expect_identical(refit_initial_equilibrium(input, current_N1), input)
  expect_identical(refit_initial_equilibrium(input, NULL), input)
  error <- tryCatch(
    refit_initial_equilibrium(input, current_N1 * 2),
    error = identity
  )
  expect_s3_class(error, "error")
  expect_match(
    conditionMessage(error),
    "Cannot safely reproduce stored initial equilibrium"
  )
})
