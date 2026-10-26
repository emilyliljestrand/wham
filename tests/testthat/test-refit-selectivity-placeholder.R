test_that("fixed legacy selectivity placeholders migrate to current arrays", {
  for (layout in c("scalar", "matrix")) {
    input <- list(
      data = list(
        n_years_model = 3L,
        n_selblocks = 2L,
        n_ages = 8L,
        selblock_models = c(2L, 1L),
        selblock_models_re = c(1L, 1L),
        selpars_est = matrix(0L, 2L, 14L),
        do_proj = 0L,
        proj_Fcatch = matrix(0)
      ),
      par = list(selpars_re = if (layout == "scalar") 0 else matrix(0)),
      map = list(selpars_re = factor(NA))
    )
    fit <- list(
      input = input,
      opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
    )
    refit <- refit_wham
    environment(refit) <- list2env(
      list(
        set_proj = identity,
        set_brp_static = function(input, ...) input,
        fit_wham = function(input, ...) list(input = input, opt = fit$opt)
      ),
      parent = environment(refit_wham)
    )
    capture.output(result <- refit(fit))
    expect_identical(result$input$par$selpars_re, array(0, c(2L, 3L, 8L)))
    expect_equal(all(is.na(result$input$map$selpars_re)), TRUE)
    expect_equal(result$input$data$selblock_models_re, c(0L, 0L))
    expect_equal(result$input$data$selpars_re_index, matrix(0L, 2L, 8L))
    expect_identical(fit$input, input)
  }
})

test_that("stored BSB Run4 fixed selectivity placeholder migrates without fitting", {
  path <- test_path("..", "..", "oldfitobjects", "BSB_2023Run4fit.rds")
  skip_if_not(file.exists(path), "Local BSB Run4 fit is unavailable")
  fit <- readRDS(path)
  refit <- refit_wham
  environment(refit) <- list2env(
    list(
      fit_wham = function(input, ...) {
        validate_input_safety(input)
        list(input = input, opt = fit$opt)
      }
    ),
    parent = environment(refit_wham)
  )
  capture.output(result <- suppressMessages(refit(fit)))
  expect_equal(dim(result$input$par$selpars_re), c(7L, 33L, 8L))
  expect_equal(all(is.na(result$input$map$selpars_re)), TRUE)
  expect_equal(result$input$data$selblock_models_re, rep(0L, 7L))
  expect_identical(result$input$par$logit_selpars, fit$parList$logit_selpars)
  expect_identical(result$input$map$logit_selpars, fit$input$map$logit_selpars)
  expect_identical(dim(result$input$par$log_NAA_sigma), c(1L, 1L, 8L))
  expect_identical(dim(result$input$par$trans_NAA_rho), c(1L, 1L, 3L))
  expect_identical(
    as.vector(result$input$par$log_NAA_sigma),
    as.vector(fit$parList$log_NAA_sigma)
  )
  expect_identical(
    as.vector(result$input$par$trans_NAA_rho),
    c(as.vector(fit$parList$trans_NAA_rho), 0)
  )
  expect_identical(result$input$map$log_NAA_sigma, fit$input$map$log_NAA_sigma)
  expect_identical(result$input$map$trans_NAA_rho, factor(rep(NA, 3L)))
  expect_identical(result$opt, fit$opt)
  for (field in c("log_NAA_sigma", "trans_NAA_rho")) {
    invalid <- result$input
    invalid$par[[field]] <- matrix(invalid$par[[field]], nrow = 1L)
    error <- tryCatch(validate_input_safety(invalid), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      paste0(field, " must have dimensions"),
      fixed = TRUE
    )
  }
})

test_that("inactive NAA matrices migrate independently of selectivity", {
  input <- list(
    data = list(
      n_stocks = 1L,
      n_regions = 1L,
      n_ages = 8L,
      n_years_model = 3L,
      NAA_re_model = 0L,
      selpars_re_index = matrix(0L)
    ),
    par = list(
      log_NAA_sigma = matrix(seq_len(8L) / 10, 1L),
      trans_NAA_rho = matrix(c(0.2, -0.3), 1L)
    ),
    map = list(
      log_NAA_sigma = factor(rep(NA, 8L), levels = "unused"),
      trans_NAA_rho = factor(rep(NA, 2L), levels = "unused")
    )
  )
  fit <- list(
    input = input,
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
  )
  fit$parList <- input$par
  fit$parList$log_NAA_sigma[] <- seq_len(8L) / 5
  refit <- refit_wham
  environment(refit) <- list2env(
    list(fit_wham = function(input, ...) list(input = input, opt = fit$opt)),
    parent = environment(refit_wham)
  )
  capture.output(result <- refit(fit))
  expect_identical(
    result$input$par$log_NAA_sigma,
    array(seq_len(8L) / 5, c(1L, 1L, 8L))
  )
  expect_identical(
    result$input$par$trans_NAA_rho,
    array(c(0.2, -0.3, 0), c(1L, 1L, 3L))
  )
  expect_identical(result$input$map$log_NAA_sigma, input$map$log_NAA_sigma)
  expect_identical(
    result$input$map$trans_NAA_rho,
    factor(rep(NA, 3L), levels = "unused")
  )
  expect_identical(fit$input, input)

  for (active in c("model", "mapping", "missing mapping", "regions")) {
    invalid <- fit
    if (active == "model") {
      invalid$input$data$NAA_re_model <- 1L
    }
    if (active == "mapping") {
      invalid$input$map$trans_NAA_rho <- factor(c(1L, NA))
    }
    if (active == "missing mapping") {
      invalid$input$map$trans_NAA_rho <- NULL
    }
    if (active == "regions") {
      invalid$input$data$n_regions <- 2L
    }
    error <- tryCatch(refit(invalid), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "Cannot safely migrate active or unsupported stored NAA matrices",
      fixed = TRUE
    )
  }
})
