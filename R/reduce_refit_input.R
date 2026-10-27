# Reject invalid peel counts before selecting terminal model years.
validate_peel_count <- function(peel, years) {
  if (
    length(peel) != 1L ||
      !is.numeric(peel) ||
      !is.finite(peel) ||
      peel != floor(peel) ||
      peel < 1L ||
      peel > length(years) - 2L
  ) {
    stop("peel must be a positive integer retaining at least two model years.")
  }
}

# Require distinct terminal years, no projections, and at least two retained years.
validate_reduction_years <- function(input, years_peeled) {
  if (input$data$n_years_proj > 0) {
    stop("cannot run wham:::reduce_input on an input for a projection")
  }
  if (
    !length(years_peeled) ||
      anyNA(years_peeled) ||
      anyDuplicated(years_peeled) ||
      !all(years_peeled %in% input$years)
  ) {
    stop("years_peeled must identify distinct model years.")
  }
  if (
    !setequal(years_peeled, tail(input$years, length(years_peeled))) ||
      length(input$years) - length(years_peeled) < 2L
  ) {
    stop(
      "reduce_input requires terminal-year peels retaining at least two model years."
    )
  }
}

# Finish a year-subset input, preserving F maps and dropping unused random fields.
reduce_refit_input <- function(new, input, ind) {
  if (!is.null(new$map$F_pars)) {
    mapping <- matrix(as.character(new$map$F_pars), nrow = NROW(new$par$F_pars))
    new$map$F_pars <- factor(mapping[ind, , drop = FALSE])
  }
  new$par$F_pars <- new$par$F_pars[ind, , drop = FALSE]
  new <- disable_empty_selblocks(new)
  new <- reduce_ecov_input(new, input)
  new$random <- Filter(
    function(field) {
      is.null(new$map[[field]]) || any(!is.na(new$map[[field]]))
    },
    new$random
  )
  new
}

# Fix all parameters and disable RE models for selectivity blocks with no years.
disable_empty_selblocks <- function(input) {
  empty_blocks <- which(input$data$n_years_selblocks == 0L)
  if (length(empty_blocks)) {
    input$data$selblock_models_re[empty_blocks] <- 0L
    for (field in c("logit_selpars", "sel_repars", "selpars_re")) {
      mapping <- input$map[[field]]
      if (is.null(mapping)) {
        mapping <- factor(seq_along(input$par[[field]]))
      }
      mapping <- array(as.character(mapping), dim(input$par[[field]]))
      if (length(dim(mapping)) == 2L) {
        mapping[empty_blocks, ] <- NA
      } else {
        mapping[empty_blocks, , ] <- NA
      }
      input$map[[field]] <- factor(mapping)
    }
  }
  input
}

# Subset Ecov years, fitted values, maps, and stored options without rebuilding.
reduce_ecov_input <- function(new, input) {
  if (!any(input$data$Ecov_model > 0) && is.null(input$options$ecov)) {
    return(new)
  }
  years <- input$years_Ecov
  if (
    !is.numeric(years) ||
      length(years) != input$data$n_years_Ecov ||
      any(!is.finite(years)) ||
      any(diff(years) != 1)
  ) {
    stop(
      "Cannot safely peel environmental covariates: years_Ecov must match the consecutive Ecov data years."
    )
  }
  keep <- which(years <= tail(new$years, 1L))
  if (!length(keep)) {
    stop("No environmental covariate years remain in this peel.")
  }
  new$years_Ecov <- years[keep]
  new$data$n_years_Ecov <- length(keep)
  new$data$years_use_Ecov <- input$data$years_use_Ecov[
    input$data$years_use_Ecov %in% (keep - 1L)
  ]
  for (field in c("Ecov_obs", "Ecov_use_obs")) {
    new$data[[field]] <- input$data[[field]][keep, , drop = FALSE]
  }
  for (field in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
    # Old inputs can have unused placeholder sigma RE matrices.
    if (NROW(input$par[[field]]) != length(years)) {
      if (
        field == "Ecov_obs_logsigma_re" &&
          !any(input$data$Ecov_obs_sigma_opt == 4L)
      ) {
        next
      }
      stop(
        "Cannot safely peel environmental covariates: invalid ",
        field,
        " dimensions."
      )
    }
    new$par[[field]] <- input$par[[field]][keep, , drop = FALSE]
    if (!is.null(input$map[[field]])) {
      mapping <- matrix(as.character(input$map[[field]]), nrow = length(years))
      new$map[[field]] <- factor(mapping[keep, , drop = FALSE])
    }
  }
  for (effect in c("R", "M", "q", "mu")) {
    field <- paste0("ind_Ecov_out_end_", effect)
    new$data[[field]] <- input$data[[field]] -
      (length(input$years) - length(new$years))
  }
  new$data$avg_years_Ecov <- input$data$avg_years_Ecov[
    input$data$avg_years_Ecov %in% (keep - 1L)
  ]
  ecov <- input$options$ecov
  if (!is.null(ecov)) {
    rows <- which(ecov$year <= tail(new$years, 1L))
    old_n <- length(ecov$year)
    ecov$year <- ecov$year[rows]
    ecov$mean <- as.matrix(ecov$mean)[rows, , drop = FALSE]
    ecov$use_obs <- as.matrix(ecov$use_obs)[rows, , drop = FALSE]
    subset_sigma <- function(x) {
      if (is.matrix(x) && nrow(x) == old_n) {
        return(x[rows, , drop = FALSE])
      }
      if (is.numeric(x) && input$data$n_Ecov == 1L && length(x) == old_n) {
        return(x[rows])
      }
      x
    }
    if (is.list(ecov$logsigma)) {
      ecov$logsigma[[1L]] <- subset_sigma(ecov$logsigma[[1L]])
    } else {
      ecov$logsigma <- subset_sigma(ecov$logsigma)
    }
    new$options$ecov <- ecov
  }
  new
}
