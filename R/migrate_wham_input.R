# Recover historical process scales only when stored likelihood components
# identify a reproducible transform. Return input plus comparison-scale flags;
# reject ambiguous fixed/shared AR1 scales rather than change the model.
refit_process_scales <- function(input, stored_report) {
  if (is.null(stored_report$nll_sel) && is.null(stored_report$nll_Ecov)) {
    return(list(input = input, selectivity = FALSE, ecov = integer()))
  }
  report_at_parameters <- function(input) {
    validate_input_safety(input)
    model <- TMB::MakeADFun(
      input$data,
      input$par,
      map = input$map,
      DLL = "wham",
      silent = TRUE
    )
    on.exit(TMB::FreeADFun(model))
    model$fn(model$par)
    model$report()
  }
  matches <- function(x, y) {
    if (is.null(x)) {
      return(all(is.finite(y) & y == 0))
    }
    isTRUE(all.equal(x, y, tolerance = 1e-8))
  }
  current <- report_at_parameters(input)
  selectivity <- FALSE
  ecov <- integer()
  if (
    !is.null(stored_report$nll_sel) &&
      !matches(current$nll_sel, stored_report$nll_sel)
  ) {
    candidate <- input
    candidate$par$sel_repars[, 2:3] <- candidate$par$sel_repars[, 2:3] * 2
    if (
      !matches(report_at_parameters(candidate)$nll_sel, stored_report$nll_sel)
    ) {
      stop("Cannot safely reproduce stored selectivity process likelihood.")
    }
    input <- candidate
    selectivity <- TRUE
  }
  if (
    !is.null(stored_report$nll_Ecov) &&
      !matches(current$nll_Ecov, stored_report$nll_Ecov)
  ) {
    ecov <- which(input$data$Ecov_model == 2L & input$data$Ecov_use_re == 1L)
    candidate <- input
    rho <- tanh(input$par$Ecov_process_pars[3, ecov] / 2)
    candidate$par$Ecov_process_pars[2, ecov] <- input$par$Ecov_process_pars[
      2,
      ecov
    ] +
      0.5 * log1p(-rho^2)
    if (
      !length(ecov) ||
        !matches(
          report_at_parameters(candidate)$nll_Ecov,
          stored_report$nll_Ecov
        )
    ) {
      stop("Cannot safely reproduce stored environmental process likelihood.")
    }
    mapping <- input$map$Ecov_process_pars
    if (!is.null(mapping)) {
      mapping <- matrix(
        as.character(mapping),
        nrow(input$par$Ecov_process_pars)
      )
      sigma <- mapping[2, ecov]
      varying_rho <- !is.na(mapping[3, ecov])
      if (
        any(
          varying_rho &
            (is.na(sigma) |
              duplicated(sigma) |
              duplicated(sigma, fromLast = TRUE))
        )
      ) {
        stop(
          "Cannot safely migrate fixed or shared historical environmental AR1 scales."
        )
      }
    }
    input <- candidate
  }
  list(input = input, selectivity = selectivity, ecov = ecov)
}

# Preserve initial-equilibrium likelihoods at stored latent values. Select the
# historical convention only when its abundance report matches the saved fit;
# stop if neither convention reproduces stored stock-level abundances.
refit_initial_equilibrium <- function(input, stored_N1) {
  stocks <- which(input$data$N1_model == 1L)
  if (!length(stocks) || is.null(stored_N1)) {
    return(input)
  }
  validate_input_safety(input)
  initial_numbers <- function(input) {
    # Compare reports at stored latent values, not newly optimized modes.
    model <- TMB::MakeADFun(
      input$data,
      input$par,
      map = input$map,
      DLL = "wham",
      silent = TRUE
    )
    on.exit(TMB::FreeADFun(model))
    model$fn(model$par)
    model$report()$N1
  }
  current <- initial_numbers(input)
  historical <- input
  historical$data$N1_model[stocks] <- -1L
  old <- initial_numbers(historical)
  for (stock in stocks) {
    target <- stored_N1[stock, , , drop = FALSE]
    if (
      isTRUE(all.equal(
        current[stock, , , drop = FALSE],
        target,
        tolerance = 1e-8
      ))
    ) {
      next
    }
    if (
      !isTRUE(all.equal(old[stock, , , drop = FALSE], target, tolerance = 1e-8))
    ) {
      stop(
        "Cannot safely reproduce stored initial equilibrium abundance for stock ",
        stock,
        "."
      )
    }
    # Preserve the historical likelihood, not just its starting abundance.
    input$data$N1_model[stock] <- -1L
  }
  input
}

# Recover recruitment decoupling from explicit options or flags first, then
# active correlation mappings. Default to coupled when evidence is absent so
# migration does not introduce an unsupported recruitment correlation.
refit_decouple_recruitment <- function(fit) {
  input <- fit$input
  options <- input$options
  NAA_re <- options[["NAA_re"]]
  if (is.null(NAA_re)) {
    NAA_re <- options[["NAA"]]
  }
  flags <- list(
    NAA_re$decouple_recruitment,
    input$data[["decouple_recruitment", exact = TRUE]],
    input$data[["decouple_rec", exact = TRUE]]
  )
  for (flag in flags) {
    if (length(flag) == 1L && !is.na(flag)) return(isTRUE(as.logical(flag)))
  }

  mapping <- input$map$trans_NAA_rho
  if (!is.null(mapping)) {
    dimensions <- dim(input$par$trans_NAA_rho)
    if (is.null(dimensions)) {
      dimensions <- c(1L, 1L, length(mapping))
    }
    if (tail(dimensions, 1L) == 3L) {
      active <- array(!is.na(mapping), dim = c(prod(head(dimensions, -1L)), 3L))
      return(any(active[, 3L]))
    }
    return(FALSE)
  }
  n_stocks <- input$data$n_stocks
  if (is.null(n_stocks)) {
    n_stocks <- 1L
  }
  n_stocks == 1L && sum(names(fit$opt$par) == "trans_NAA_rho") == 3L
}

# Migrate supported single-stock pre-2.0 fits, retaining fitted values, maps,
# observations, and correlation conventions. Reject unsupported structures or
# dimension mismatches instead of silently rebuilding a different model.
migrate_legacy_refit_input <- function(fit) {
  old <- fit$input
  data <- old$data
  data$decouple_rec <- as.integer(refit_decouple_recruitment(fit))
  unsupported <- character()
  if (is.null(old$asap3$n_ages)) {
    unsupported <- c(unsupported, "missing single-stock ASAP3 data")
  }
  if (any(data$Ecov_model != 0)) {
    unsupported <- c(unsupported, "environmental covariates")
  }
  if (any(data$selblock_models_re != 1)) {
    unsupported <- c(unsupported, "selectivity random effects")
  }
  if (any(data$M_re_model != 1)) {
    unsupported <- c(unsupported, "mortality random effects")
  }
  if (data$M_model == 3) {
    unsupported <- c(unsupported, "weight-based mortality")
  }
  if (isTRUE(data$use_steepness == 1)) {
    unsupported <- c(unsupported, "steepness parameterization")
  }
  if (isTRUE(data$do_proj == 1)) {
    unsupported <- c(unsupported, "stored projections")
  }
  if (length(unsupported)) {
    stop(
      "Cannot safely migrate this legacy WHAM input: ",
      paste(unsupported, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  asap3 <- align_fitted_asap3_indices(old$asap3)
  input <- prepare_wham_input(
    asap3 = list(list(dat = asap3)),
    model_name = old$model_name,
    recruit_model = data$recruit_model,
    NAA_re = list(decouple_recruitment = as.logical(data$decouple_rec))
  )
  for (field in intersect(names(data), names(input$data))) {
    value <- data[[field]]
    target <- input$data[[field]]
    if (is.null(dim(target)) || is.data.frame(target) || field == "waa") {
      input$data[[field]] <- value
    } else if (length(value) == length(target)) {
      target[] <- value
      input$data[[field]] <- target
    }
  }
  input$data$waa_pointer_indices <- data$waa_pointer_indices[seq_len(
    data$n_indices
  )]
  input$data$waa_pointer_M <- data$waa_pointer_jan1
  input$data$M_model <- 1L
  input$data$selblock_models_re <- rep(0L, data$n_selblocks)
  input$data$NAA_re_model <- min(data$n_NAA_sigma, 2L)
  input$data$decouple_recruitment <- data$decouple_rec
  input$data$use_alt_AR1 <- 0L
  if (!is.null(data$static_FXSPR_init)) {
    input$data$FXSPR_static_init <- data$static_FXSPR_init
  }
  if (!is.null(data$avg_years_ind)) {
    for (field in grep("^avg_years_ind_", names(input$data), value = TRUE)) {
      target <- input$data[[field]]
      if (length(dim(target)) == 2) {
        target[,] <- 0
        target[seq_len(length(data$avg_years_ind) + 1L), ] <- c(
          length(data$avg_years_ind),
          data$avg_years_ind
        )
      } else if (length(dim(target)) == 3) {
        target[] <- 0
        target[,, seq_len(length(data$avg_years_ind) + 1L)] <- rep(
          c(length(data$avg_years_ind), data$avg_years_ind),
          each = prod(dim(target)[1:2])
        )
      }
      input$data[[field]] <- target
    }
  }

  pars <- old$par
  if (!is.null(fit$parList)) {
    pars[names(fit$parList)] <- fit$parList
  }
  transfer <- function(
    target,
    source = target,
    indices = seq_along(pars[[source]]),
    scale = 1
  ) {
    values <- pars[[source]][indices] * scale
    if (length(values) != length(input$par[[target]])) {
      stop(
        "Cannot safely migrate legacy parameter dimensions: ",
        source,
        call. = FALSE
      )
    }
    input$par[[target]][] <<- values
    mapping <- old$map[[source]]
    if (is.null(mapping)) {
      mapping <- factor(seq_along(pars[[source]]))
    }
    input$map[target] <<- list(factor(
      as.character(mapping[indices]),
      levels = levels(mapping)
    ))
  }
  renamed <- c(
    "log_N1_pars",
    "log_NAA_sigma",
    "mean_rec_pars",
    "F_devs",
    "log_F1",
    "M_a"
  )
  inactive <- c("selpars_re", "Ecov_beta", "Ecov_obs_sigma_par")
  for (parameter in setdiff(names(pars), c(renamed, inactive))) {
    if (parameter %in% names(input$par)) transfer(parameter)
  }
  for (parameter in c("log_N1", "mean_rec_pars")) {
    source <- if (parameter == "log_N1") "log_N1_pars" else parameter
    values <- pars[[source]]
    input$par[[parameter]][seq_along(values)] <- values
    mapping <- rep(NA_integer_, length(input$par[[parameter]]))
    source_map <- old$map[[source]]
    if (is.null(source_map)) {
      source_map <- factor(seq_along(values))
    }
    mapping[seq_along(values)] <- as.integer(source_map)
    input$map[[parameter]] <- factor(mapping)
  }
  pointers <- data$NAA_sigma_pointers
  if (data$n_NAA_sigma > 0) {
    transfer("log_NAA_sigma", indices = pointers)
  }
  transfer(
    "Mpars",
    "M_a",
    if (data$M_model == 1) rep(1L, data$n_ages) else seq_len(data$n_ages)
  )
  input$par$trans_NAA_rho[] <- 2 * pars$trans_NAA_rho
  if (data$decouple_rec == 1) {
    input$par$trans_NAA_rho[] <- 2 * pars$trans_NAA_rho[c(1L, 3L, 2L)]
    mapping <- input$map$trans_NAA_rho
    input$map$trans_NAA_rho <- factor(
      as.character(mapping[c(1L, 3L, 2L)]),
      levels = levels(mapping)
    )
    if (data$n_NAA_sigma > 1) {
      rho <- tanh(pars$trans_NAA_rho)
      input$par$log_NAA_sigma[1] <- input$par$log_NAA_sigma[1] -
        0.5 * log1p(-rho[1]^2)
      input$par$log_NAA_sigma[-1] <- input$par$log_NAA_sigma[-1] +
        0.5 * log1p(-rho[3]^2) -
        0.5 * log1p(-rho[2]^2)
    }
  }
  input$par$q_repars[, 2] <- 2 * pars$q_repars[, 2]
  input$par$sel_repars[, 2:3] <- 2 * pars$sel_repars[, 2:3]
  input$par$M_repars[1, 1, 2:3] <- 2 * pars$M_repars[2:3]
  input$par$F_pars[1, ] <- pars$log_F1
  input$par$F_pars[-1, ] <- pars$F_devs
  map_F <- matrix(NA_character_, data$n_years_model, data$n_fleets)
  for (parameter in c("log_F1", "F_devs")) {
    mapping <- old$map[[parameter]]
    if (is.null(mapping)) {
      mapping <- factor(seq_along(pars[[parameter]]))
    }
    values <- ifelse(
      is.na(mapping),
      NA_character_,
      paste(parameter, as.character(mapping))
    )
    if (parameter == "log_F1") map_F[1, ] <- values else map_F[-1, ] <- values
  }
  input$map$F_pars <- factor(
    map_F,
    levels = unique(as.vector(map_F)[!is.na(map_F)])
  )
  input$random <- old$random
  input$options <- NULL
  input
}
