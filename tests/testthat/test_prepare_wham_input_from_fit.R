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
  input$data$decouple_recruitment <- NULL

  rebuilt <- prepare_wham_input_from_fit(list(input = input))

  expect_identical(rebuilt$data$decouple_recruitment, 0L)
})

test_that("invalid fitted model files produce an informative error", {
  expect_error(
    prepare_wham_input_from_fit(list()),
    "does not contain an \\$input component"
  )
})

test_that("recruitment decoupling uses stored flags and active correlations", {
  fit <- list(input = list(data = list(), par = list(trans_NAA_rho = numeric(3))))
  expect_false(refit_decouple_recruitment(fit))
  fit$opt$par <- setNames(numeric(3), rep("trans_NAA_rho", 3))
  expect_true(refit_decouple_recruitment(fit))
  fit$input$map$trans_NAA_rho <- factor(c(1L, 2L, NA_integer_))
  expect_false(refit_decouple_recruitment(fit))
  fit$input$map$trans_NAA_rho <- factor(1:3)
  expect_true(refit_decouple_recruitment(fit))
  fit$input$data$decouple_rec <- 0L
  expect_false(refit_decouple_recruitment(fit))
  fit$input$data$decouple_rec <- 1L
  expect_true(refit_decouple_recruitment(fit))
  fit$input$options$NAA_re$decouple_recruitment <- FALSE
  expect_false(refit_decouple_recruitment(fit))
})

test_that("rebuilding preserves decoupling when preparation options omit it", {
  asap3 <- read_asap3_dat(
    file.path(system.file("extdata", package = "wham"), "ex1_SNEMAYT.dat")
  )
  input <- prepare_wham_input(asap3)
  input$options$NAA_re$decouple_recruitment <- NULL
  input$data$decouple_recruitment <- 1L
  expect_identical(prepare_wham_input_from_fit(list(input = input))$data$decouple_recruitment, 1L)
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

  capture.output(result <- refit(fit))
  expect_identical(result$input$data$years_use, 0:2)
  expect_identical(result$input$data$decouple_recruitment, 0L)
  expect_null(fit$input$data$years_use)

  fit$input$data$years_use <- c(0L, 2L)
  fit$input$data$decouple_recruitment <- 1L
  capture.output(result <- refit(fit))
  expect_identical(result$input, fit$input)
})

test_that("refitting preserves stored mappings and uses fitted starting values", {
  fit <- list(
    input = list(
      data = list(n_years_model = 3L),
      par = list(parameter = c(0, 7, 0)),
      map = list(parameter = factor(c(1, NA, 1))),
      random = character(),
      options = list()
    ),
    parList = list(parameter = c(2, 7, 2)),
    opt = list(par = c(parameter = 2), objective = 0, convergence = 0)
  )
  original_input <- fit$input
  refit <- refit_wham
  environment(refit) <- list2env(list(
    prepare_wham_input_from_fit = function(...) stop("Input must not be rebuilt"),
    fit_wham = function(input, ...) {
      result <- fit
      result$input <- input
      result
    }
  ), parent = environment(refit_wham))
  capture.output(result <- refit(fit))
  expect_identical(result$input$map, original_input$map)
  expect_identical(result$input$random, original_input$random)
  expect_identical(result$input$par, fit$parList)
  expect_identical(fit$input, original_input)
})

test_that("packed and array selectivity migrations retain means and shared effects", {
  for (layout in c("packed", "array")) for (model in c(4L, 5L)) {
    shared <- model == 4L
    values <- if (shared) c(0.1, 0.2, 0.1, 0.2) else seq_len(4) / 10
    mapping <- if (shared) c(2, 1, 2, 1) else c(4, 3, 2, 1)
    deviations <- values
    if (layout == "array") {
      deviations <- array(0, c(1L, 2L, 4L))
      deviations[1, , 1:2] <- values
      array_mapping <- array(NA_integer_, dim(deviations))
      array_mapping[1, , 1:2] <- mapping
      mapping <- array_mapping
    }
    input <- list(
      data = list(
        n_years_model = 2L, n_selblocks = 1L, n_ages = 2L,
        selblock_models = 1L, selblock_models_re = model,
        selblock_years = matrix(1L, 2L, 1L), n_selpars_est = 2L,
        selpars_est = matrix(c(1L, 1L, rep(0L, 6)), 1L),
        selpars_lower = matrix(0, 1L, 8L), selpars_upper = matrix(1, 1L, 8L),
        do_proj = 0L, proj_Fcatch = matrix(0)
      ),
      par = list(logit_selpars = matrix(seq_len(8), 1L), selpars_re = deviations),
      map = list(logit_selpars = factor(c(2, 1, rep(NA, 6))), selpars_re = factor(mapping))
    )
    fit <- list(input = input, opt = list(par = c(parameter = 1), objective = 0, convergence = 0))
    refit <- refit_wham
    environment(refit) <- list2env(list(
      set_proj = identity,
      set_brp_static = function(input, ...) input,
      fit_wham = function(input, ...) list(input = input, opt = fit$opt)
    ), parent = environment(refit_wham))
    capture.output(result <- refit(fit))
    migrated <- result$input
    expect_identical(migrated$par$logit_selpars, input$par$logit_selpars)
    expect_identical(migrated$map$logit_selpars, input$map$logit_selpars)
    expect_identical(migrated$data$selpars_lower, input$data$selpars_lower)
    expect_identical(migrated$data$selpars_upper, input$data$selpars_upper)
    expect_equal(migrated$data$n_selpars_re, if (shared) 1L else 2L)
    expect_equal(migrated$data$selpars_re_index[1, 1:2], if (shared) c(1L, 1L) else 1:2)
    expect_equal(as.vector(migrated$par$selpars_re[1, , seq_len(if (shared) 1L else 2L)]),
                 if (shared) values[1:2] else values)
    expect_identical(levels(migrated$map$selpars_re), levels(input$map$selpars_re))
  }
})

test_that("stored BSB selectivity arrays preserve parameters and likelihood", {
  fit_file <- test_path("..", "..", "oldfitobjects", "BSB_2024MTfit.RDS")
  skip_if_not(file.exists(fit_file), "Local BSB fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) {
      expect_identical(input$par$logit_selpars, fit$parList$logit_selpars)
      expect_identical(input$map$logit_selpars, fit$input$map$logit_selpars)
      expect_identical(input$map$log_index_sig_scale, fit$input$map$log_index_sig_scale)
      expect_identical(input$data$obsvec, fit$input$data$obsvec)
      model <- TMB::MakeADFun(
        input$data, input$par, DLL = "wham", random = input$random,
        map = input$map, silent = TRUE
      )
      expect_equal(model$par, fit$opt$par)
      expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)
      expect_lt(max(abs(model$gr(model$par))), 1e-6)
      report <- model$report()
      for (field in c("SSB", "NAA", "pred_catch", "pred_indices", "nll_sel")) {
        expect_equal(as.vector(report[[field]]), as.vector(fit$rep[[field]]), tolerance = 1e-8)
      }
      list(input = input, opt = fit$opt, parList = input$par)
    }
  ), parent = environment(refit_wham))
  capture.output(result <- suppressMessages(refit(fit)))
  expect_length(result$opt$par, 211L)
  expect_identical(fit$input, original_input)
})

test_that("legacy migration defaults missing recruitment information before rebuilding", {
  input <- list(
    asap3 = list(n_ages = 5L),
    data = list(
      Ecov_model = 0, selblock_models_re = 1, M_re_model = 1,
      M_model = 2, use_steepness = 0, do_proj = 0,
      decouple_rec_unused = 1L
    )
  )
  migrate <- migrate_legacy_refit_input
  environment(migrate) <- list2env(list(
    prepare_wham_input = function(..., NAA_re) {
      expect_identical(NAA_re$decouple_recruitment, FALSE)
      stop("Reached input rebuilding")
    }
  ), parent = environment(migrate_legacy_refit_input))

  expect_error(migrate(list(input = input)), "Reached input rebuilding")
  expect_null(input$data[["decouple_rec", exact = TRUE]])
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
    M_model = 3, use_steepness = 1, do_proj = 1
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

test_that("legacy decoupled Cod migration preserves the likelihood", {
  fit_file <- test_path("..", "..", "oldfitobjects", "GBK_Cod_2024MTfit.rds")
  skip_if_not(file.exists(fit_file), "Local decoupled fit is unavailable")
  fit <- readRDS(fit_file)
  original_input <- fit$input
  input <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(fit$input, original_input)
  expect_identical(input$data$decouple_recruitment, 1L)
  expect_equal(as.vector(input$par$trans_NAA_rho), 2 * fit$parList$trans_NAA_rho[c(1L, 3L, 2L)])

  model <- TMB::MakeADFun(
    input$data, input$par, DLL = "wham", random = input$random,
    map = input$map, silent = TRUE
  )
  expect_equal(as.numeric(model$fn(model$par)), as.numeric(fit$opt$objective), tolerance = 1e-8)

  refit <- refit_wham
  environment(refit) <- list2env(list(
    fit_wham = function(input, ...) list(input = input, opt = fit$opt, parList = input$par)
  ), parent = environment(refit_wham))
  expect_message(
    output <- capture.output(result <- suppressMessages(refit(fit))),
    NA
  )
  expect_identical(result$input$data$decouple_recruitment, 1L)
  expect_true(any(grepl("max_abs_par_diff", output)))

  fit$input$data$decouple_rec <- NULL
  inferred <- suppressMessages(migrate_legacy_refit_input(fit))
  expect_identical(inferred$data$decouple_recruitment, 1L)
  expect_equal(inferred$par, input$par)
  expect_identical(inferred$map, input$map)
})
