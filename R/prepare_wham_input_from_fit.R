#' Prepare input data and parameters from a fitted WHAM model
#'
#' Rebuilds input using the preparation options stored in a fitted WHAM model.
#' The current \code{\link{prepare_wham_input}} applies known input-structure
#' changes to older fitted models.
#'
#' @param fit a fitted WHAM model object containing an \code{$input} component
#'
#' @return A named list suitable for use as the \code{input} argument to
#'   \code{\link{fit_wham}}. For fitted models created before preparation
#'   options were stored, the original input is returned unchanged with a
#'   warning.
#'
#' @seealso \code{\link{prepare_wham_input}}, \code{\link{fit_wham}}
#'
#' @export
prepare_wham_input_from_fit <- function(fit = NULL) {
  if (is.null(fit)) {
    stop("Provide a fitted WHAM model object.")
  }
  if (!is.list(fit) || is.null(fit$input)) {
    stop("The fitted WHAM model does not contain an $input component.")
  }
  old_input <- fit$input
  if (!is.list(old_input)) {
    stop("The fitted model's $input component is not a list.")
  }

  options <- old_input$options
  if (is.null(options) || !is.list(options)) {
    warning(
      "The fitted model does not contain preparation options; returning its stored input unchanged."
    )
    return(old_input)
  }

  # Translate historical option names (e.g. "NAA" vs "NAA_re") here, taking the
  # first candidate present without changing names used by component setters.
  option_names <- names(options)
  pick <- function(...) {
    hits <- intersect(c(...), option_names)
    if (length(hits)) options[[hits[1]]] else NULL
  }
  args <- list(
    basic_info = pick("basic_info"),
    catch_info = pick("catch_info", "catch"),
    index_info = pick("index_info", "index"),
    waa = pick("waa"),
    NAA_re = pick("NAA_re", "NAA"),
    catchability = pick("catchability", "q"),
    selectivity = pick("selectivity"),
    age_comp = pick("age_comp"),
    F = pick("F"),
    M = pick("M"),
    move = pick("move"),
    L = pick("L"),
    ecov = pick("ecov")
  )

  if (is.null(args$NAA_re$decouple_recruitment)) {
    if (is.null(args$NAA_re)) {
      args$NAA_re <- list()
    }
    args$NAA_re$decouple_recruitment <- refit_decouple_recruitment(fit)
  }

  # WAA options were historically stored separately from basic_info, while
  # prepare_wham_input receives them through basic_info.
  if (!is.null(args$waa)) {
    if (is.null(args$basic_info)) {
      args$basic_info <- list()
    }
    args$basic_info[names(args$waa)] <- args$waa
  }
  args$waa <- NULL

  if (is.null(args$basic_info)) {
    args$basic_info <- list()
  }
  recruit_model <- args$basic_info$recruit_model
  if (is.null(recruit_model)) {
    recruit_model <- old_input$data$recruit_model
  }
  if (is.null(recruit_model)) {
    recruit_model <- 2
  }

  # prepare_wham_input expects the structure returned by read_asap3_dat.
  # Fitted inputs retain the processed $dat elements, so restore that wrapper.
  asap3 <- old_input$asap3
  if (!is.null(asap3)) {
    if (!is.list(asap3)) {
      stop("The fitted model's $asap3 component is not a list.")
    }
    if (
      !all(vapply(asap3, function(x) is.list(x) && !is.null(x$dat), logical(1)))
    ) {
      asap3 <- lapply(asap3, function(x) list(dat = x))
    }
  }

  if (!is.null(asap3)) {
    asap3 <- lapply(asap3, function(x) {
      x$dat <- align_fitted_asap3_indices(x$dat)
      x
    })
  }

  args$asap3 <- asap3
  args$model_name <- old_input$model_name
  args$recruit_model <- recruit_model
  args <- args[!vapply(args, is.null, logical(1))]

  rebuilt <- tryCatch(
    do.call(prepare_wham_input, args),
    error = function(error) {
      stop(
        "Could not rebuild input from the fitted model's preparation options: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  # Preserve structural fields required by TMB when rebuilding older fits.
  if (is.null(rebuilt$data$years_use) && !is.null(old_input$data$years_use)) {
    rebuilt$data$years_use <- old_input$data$years_use
  }
  rebuilt
}

# Fitted models store ASAP3 index data after subsetting IAA_mats and the
# per-index fields to the used indices, while use_index keeps its full length.
# Realign those fields so set_indices() can select the used indices again.
align_fitted_asap3_indices <- function(dat) {
  if (is.null(dat$use_index) || is.null(dat$IAA_mats)) {
    return(dat)
  }
  n_full <- length(dat$use_index)
  n_stored <- length(dat$IAA_mats)
  if (n_stored == n_full || n_stored != sum(dat$use_index == 1)) {
    return(dat)
  }
  keep <- which(dat$use_index == 1)
  index_fields <- c(
    "index_units",
    "index_acomp_units",
    "index_month",
    "use_index_acomp",
    "index.names",
    "index_WAA_pointers",
    "index_sel_choice",
    "index_sel_option",
    "index_sel_start_age",
    "index_sel_end_age",
    "index_sel_ini"
  )
  for (f in index_fields) {
    if (length(dat[[f]]) == n_full) dat[[f]] <- dat[[f]][keep]
  }
  dat$use_index <- rep(1L, n_stored)
  dat
}
