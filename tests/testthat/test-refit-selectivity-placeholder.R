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

  for (active in c("model", "mapping", "missing mapping")) {
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
    error <- tryCatch(refit(invalid), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "Cannot safely migrate active or unsupported stored NAA matrices",
      fixed = TRUE
    )
  }
})

test_that("inactive NAA matrices expand by stock and region without fitting", {
  input <- list(
    data = list(
      n_stocks = 2L,
      n_regions = 2L,
      n_ages = 3L,
      n_years_model = 3L,
      NAA_re_model = c(0L, 0L),
      mu_model = 1,
      can_move = array(0, c(2L, 1L, 2L, 2L)),
      selpars_re_index = matrix(0L)
    ),
    par = list(
      log_NAA_sigma = matrix(seq_len(6L) / 10, 2L),
      trans_NAA_rho = matrix(c(0.1, 0.2, -0.3, -0.4), 2L)
    ),
    map = list(
      log_NAA_sigma = factor(rep(NA, 6L), levels = "unused"),
      trans_NAA_rho = factor(rep(NA, 4L), levels = "unused")
    )
  )
  fit <- list(
    input = input,
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
  )
  refit <- refit_wham
  environment(refit) <- list2env(
    list(fit_wham = function(input, ...) list(input = input, opt = fit$opt)),
    parent = environment(refit_wham)
  )
  capture.output(result <- refit(fit))
  expect_identical(result$input$data$mu_model, matrix(1L, 2L, 1L))
  expect_identical(dim(result$input$par$log_NAA_sigma), c(2L, 2L, 3L))
  expect_identical(dim(result$input$par$trans_NAA_rho), c(2L, 2L, 3L))
  for (region in 1:2) {
    expect_equal(
      result$input$par$log_NAA_sigma[, region, ],
      input$par$log_NAA_sigma
    )
    expect_equal(
      result$input$par$trans_NAA_rho[, region, 1:2],
      input$par$trans_NAA_rho
    )
  }
  expect_equal(result$input$par$trans_NAA_rho[,, 3L], matrix(0, 2L, 2L))
  expect_identical(
    result$input$map$log_NAA_sigma,
    factor(rep(NA, 12L), levels = "unused")
  )
  expect_identical(
    result$input$map$trans_NAA_rho,
    factor(rep(NA, 12L), levels = "unused")
  )
  expect_identical(fit$input, input)
  fit$input$data$can_move[1L] <- 1
  capture.output(result <- refit(fit))
  expect_identical(result$input$data$mu_model, 1)
})

test_that("independent active NAA matrices preserve fitted values and sigma sharing", {
  input <- list(
    data = list(
      n_stocks = 1L,
      n_regions = 1L,
      n_ages = 3L,
      n_years_model = 3L,
      NAA_re_model = 2L,
      selpars_re_index = matrix(0L)
    ),
    par = list(
      log_NAA_sigma = matrix(c(-0.1, -0.3, -0.3), 1L),
      trans_NAA_rho = matrix(0, 1L, 2L)
    ),
    map = list(
      log_NAA_sigma = factor(c(1L, 2L, 2L)),
      trans_NAA_rho = factor(c(NA, NA))
    ),
    random = "log_NAA"
  )
  fit <- list(
    input = input,
    parList = input$par,
    opt = list(par = c(parameter = 1), objective = 0, convergence = 0)
  )
  fit$parList$log_NAA_sigma[] <- c(-0.2, -0.4, -0.4)
  refit <- refit_wham
  environment(refit) <- list2env(
    list(fit_wham = function(input, ...) list(input = input, opt = fit$opt)),
    parent = environment(refit_wham)
  )
  capture.output(result <- refit(fit))
  expect_identical(
    result$input$par$log_NAA_sigma,
    array(c(-0.2, -0.4, -0.4), c(1L, 1L, 3L))
  )
  expect_identical(result$input$map$log_NAA_sigma, input$map$log_NAA_sigma)
  expect_identical(result$input$par$trans_NAA_rho, array(0, c(1L, 1L, 3L)))
  expect_identical(result$input$map$trans_NAA_rho, factor(rep(NA, 3L)))
  expect_identical(result$input$random, "log_NAA")
  expect_identical(fit$input, input)

  for (case in c(
    "nonzero",
    "estimated zero",
    "missing map",
    "short map"
  )) {
    invalid <- fit
    if (case == "nonzero") {
      invalid$parList$trans_NAA_rho[1L] <- 0.1
    }
    if (case == "estimated zero") {
      invalid$input$map$trans_NAA_rho <- factor(c(1L, NA))
    }
    if (case == "missing map") {
      invalid$input$map$log_NAA_sigma <- NULL
    }
    if (case == "short map") {
      invalid$input$map$log_NAA_sigma <- factor(1L)
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

test_that("active multiregion NAA migration retains shared sigma groups", {
  files <- c("BSB_2023Run11Combinedfit.rds", "BSB_2023Run17Combinedfit.rds")
  for (file in files) {
    path <- test_path("..", "..", "oldfitobjects", file)
    skip_if_not(file.exists(path), "Local combined fits are unavailable")
    fit <- readRDS(path)
    fit$call <- quote(fit_wham(
      input = input,
      do.retro = FALSE,
      do.osa = FALSE,
      do.sdrep = FALSE,
      do.brps = FALSE,
      MakeADFun.silent = TRUE
    ))
    capture.output(result <- suppressMessages(refit_wham(fit)))
    expect_lt(abs(result$opt$objective - fit$opt$objective), 1e-6)
    expect_equal(result$opt$convergence, 0)
    expect_identical(length(result$opt$par), length(fit$opt$par))
    mapping <- array(
      as.character(result$input$map$log_NAA_sigma),
      dim(result$input$par$log_NAA_sigma)
    )
    for (region in seq_len(fit$input$data$n_regions)) {
      expect_identical(
        as.vector(mapping[, region, ]),
        as.character(fit$input$map$log_NAA_sigma)
      )
      expect_equal(
        result$input$par$log_NAA_sigma[, region, ],
        fit$parList$log_NAA_sigma
      )
    }
  }
})

test_that("six stored matrix-layout fits reproduce their original likelihoods", {
  files <- c(
    "north_south_scaa_fit_less_good.rds",
    "bridge9_combinedfit.rds",
    "north_runs_fit_bestfit.rds",
    "tfit_d_d.rds",
    "tfit_dm_dm_vast_constant.rds",
    "tfit_dm_dm.rds"
  )
  for (file in files) {
    path <- test_path("..", "..", "oldfitobjects", file)
    skip_if_not(file.exists(path), "Local legacy fits are unavailable")
    fit <- readRDS(path)
    fit$call <- quote(fit_wham(
      input = input,
      do.retro = FALSE,
      do.osa = FALSE,
      do.sdrep = FALSE,
      do.brps = FALSE,
      MakeADFun.silent = TRUE
    ))
    capture.output(result <- suppressMessages(refit_wham(fit)))
    expect_lt(abs(result$opt$objective - fit$opt$objective), 1e-6)
    expect_equal(result$opt$convergence, 0)
    data <- fit$input$data
    expect_identical(
      dim(result$input$par$log_NAA_sigma),
      as.integer(c(data$n_stocks, data$n_regions, data$n_ages))
    )
    expect_identical(
      dim(result$input$par$trans_NAA_rho),
      as.integer(c(data$n_stocks, data$n_regions, 3L))
    )
    expect_equal(all(is.na(result$input$map$trans_NAA_rho)), TRUE)
    if (data$n_regions == 1L) {
      expect_identical(
        result$input$map$log_NAA_sigma,
        fit$input$map$log_NAA_sigma
      )
      expect_identical(
        as.vector(result$input$par$log_NAA_sigma),
        as.vector(fit$parList$log_NAA_sigma)
      )
    }
  }
})
