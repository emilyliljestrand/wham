validate_input_safety <- function(input) {
  data <- input$data
  par <- input$par
  fail <- function(...) stop("Unsafe WHAM input: ", ..., call. = FALSE)
  integer_values <- function(x) {
    is.numeric(x) && all(is.finite(x) & x == floor(x))
  }
  shape <- function(x, expected, field) {
    if (!identical(as.integer(dim(x)), as.integer(expected))) {
      fail(
        field,
        " must have dimensions ",
        paste(expected, collapse = " x "),
        "."
      )
    }
  }
  ny <- data$n_years_model
  ne <- data$n_Ecov
  nye <- data$n_years_Ecov
  for (field in c("n_years_model", "n_Ecov", "n_years_Ecov")) {
    x <- data[[field]]
    if (length(x) != 1L || !integer_values(x) || x < 1L) {
      fail(field, " must be a positive integer.")
    }
  }
  if (ny < 2L) {
    fail("at least two model years must remain.")
  }
  for (field in names(input$map)) {
    mapping <- input$map[[field]]
    if (
      !is.factor(mapping) ||
        is.null(par[[field]]) ||
        length(mapping) != length(par[[field]])
    ) {
      fail("map$", field, " must be a factor matching par$", field, ".")
    }
  }
  if (length(setdiff(input$random, names(par)))) {
    fail("random contains unknown parameters.")
  }
  for (field in c("Ecov_model", "Ecov_use_re", "Ecov_obs_sigma_opt")) {
    if (length(data[[field]]) != ne || !integer_values(data[[field]])) {
      fail(field, " must contain n_Ecov integers.")
    }
  }
  options <- list(
    Ecov_model = 0:2,
    Ecov_use_re = 0:1,
    Ecov_obs_sigma_opt = c(1L, 2L, 4L),
    Ecov_use_obs = 0:1
  )
  for (field in names(options)) {
    if (!all(data[[field]] %in% options[[field]])) {
      fail(field, " contains invalid option codes.")
    }
  }
  shape(par$Ecov_process_pars, c(3L, ne), "Ecov_process_pars")
  if (
    !is.matrix(par$Ecov_re) ||
      ncol(par$Ecov_re) != ne ||
      nrow(par$Ecov_re) < nye
  ) {
    fail("Ecov_re must span all environmental data years and covariates.")
  }
  for (field in c("Ecov_obs", "Ecov_use_obs")) {
    shape(data[[field]], c(nye, ne), field)
  }
  shape(par$Ecov_obs_logsigma, c(nye, ne), "Ecov_obs_logsigma")
  if (any(data$Ecov_obs_sigma_opt == 4L)) {
    shape(par$Ecov_obs_logsigma_re, c(nye, ne), "Ecov_obs_logsigma_re")
    shape(par$Ecov_obs_sigma_par, c(2L, ne), "Ecov_obs_sigma_par")
    cols <- which(data$Ecov_obs_sigma_opt == 4L)
    if (
      any(!is.finite(par$Ecov_obs_logsigma_re[, cols])) ||
        any(!is.finite(par$Ecov_obs_sigma_par[, cols]))
    ) {
      fail(
        "environmental observation-error random effects need finite starting values."
      )
    }
  }
  used <- data$Ecov_use_obs == 1L
  if (anyNA(used) || any(!is.finite(data$Ecov_obs[used]))) {
    fail(
      "used environmental observations must be finite; mark missing observations unused."
    )
  }
  fixed_sigma <- used
  fixed_sigma[, data$Ecov_obs_sigma_opt == 4L] <- FALSE
  if (any(!is.finite(par$Ecov_obs_logsigma[fixed_sigma]))) {
    fail("used environmental observations need finite log standard deviations.")
  }
  indices <- data$years_use_Ecov
  if (
    !integer_values(indices) ||
      any(indices < 0L | indices >= nrow(par$Ecov_re)) ||
      (any(data$Ecov_model > 0L) && !length(indices))
  ) {
    fail("years_use_Ecov contains invalid environmental year indices.")
  }
  for (effect in c("R", "M", "q", "mu")) {
    how <- data[[paste0("Ecov_how_", effect)]]
    if (
      !integer_values(how) ||
        any(how < 0L | how > if (effect == "R") 5L else 1L)
    ) {
      fail("Ecov_how_", effect, " contains invalid effect codes.")
    }
    if (!any(how > 0L)) {
      next
    }
    start <- data[[paste0("ind_Ecov_out_start_", effect)]]
    end <- data[[paste0("ind_Ecov_out_end_", effect)]]
    degree <- data[[paste0("n_poly_Ecov_", effect)]]
    for (field in c("start", "end", "degree")) {
      x <- get(field)
      if (!integer_values(x) || !identical(dim(x), dim(how))) {
        fail("invalid Ecov ", effect, " ", field, " dimensions or values.")
      }
    }
    # Native code evaluates every covariate in a family if any effect is active.
    if (any(start < 0L | end >= nrow(par$Ecov_re) | end - start + 1L != ny)) {
      fail(
        "Ecov ",
        effect,
        " start/end indices must be in bounds and span exactly n_years_model; rebuild or reduce the covariate input."
      )
    }
    if (any(degree < 1L | degree >= ny)) {
      fail(
        "Ecov ",
        effect,
        " polynomial degree must be positive and smaller than the retained model length."
      )
    }
  }
  nb <- data$n_selblocks
  if (length(nb) != 1L || !integer_values(nb) || nb < 1L) {
    fail("n_selblocks must be a positive integer.")
  }
  shape(data$selblock_years, c(ny, nb), "selblock_years")
  if (
    anyNA(data$selblock_years) ||
      !all(data$selblock_years %in% 0:1) ||
      !isTRUE(all.equal(
        as.numeric(colSums(data$selblock_years)),
        as.numeric(data$n_years_selblocks)
      ))
  ) {
    fail("n_years_selblocks must equal the column sums of selblock_years.")
  }
  shape(par$selpars_re, c(nb, ny, max(data$n_ages, 4L)), "selpars_re")
  shape(par$sel_repars, c(nb, 3L), "sel_repars")
  active <- data$selblock_models_re > 0L
  if (
    length(active) != nb ||
      anyNA(active) ||
      length(data$n_selpars_re) != nb ||
      !integer_values(data$n_selpars_re) ||
      any(
        data$n_selpars_re[active] < 1L |
          data$n_selpars_re[active] > dim(par$selpars_re)[3L]
      )
  ) {
    fail("invalid selectivity random-effect counts.")
  }
  if (any(active & data$n_years_selblocks == 0L)) {
    fail("active selectivity random-effect block has no retained years.")
  }
  invisible(input)
}
