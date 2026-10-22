# Manual, opt-in refit checks for saved WHAM fits.
# Change one setting at a time. The path is resolved from either the project
# root or the tests directory, so the working directory does not matter.

if (!requireNamespace("wham", quietly = TRUE)) {
  stop("Install or load the wham package before running this script.")
}

# Load the local package source, including the full repository's R code.
devtools::load_all()

# Change this setting before each run.
# fit_name <- "BSB_2025MTfit.rds.RDS"
fit_name <- "BSB_2024MTfit.rds.RDS"
# fit_name <- "BSB_2023RTfit.rds.RDS"
# fit_name <- "GBKCod_2025MTfit.rds.rds"
check_name <- c(
  "check_convergence",
  "aic",
  "mohns_rho",
  "compare_wham_models",
  "plot_wham_output",
  "self_test"
)
run_refit <- TRUE

project_dirs <- normalizePath(
  c(getwd(), file.path(getwd(), "..")),
  mustWork = FALSE
)
project_dir <- project_dirs[dir.exists(file.path(
  project_dirs,
  "oldfitobjects"
))][1]
if (is.na(project_dir)) stop("Could not locate the oldfitobjects directory.")

fit_file <- file.path(project_dir, "oldfitobjects", fit_name)
if (!file.exists(fit_file)) stop("Fit file does not exist: ", fit_file)

# Use the test session's temporary directory for generated output.
output_dir <- file.path(tempdir(), "test_refit_wham")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and inspect one saved fit object.
fit <- readRDS(fit_file)
tools::file_path_sans_ext(basename(fit_file))

# Errors are deliberately not caught so they are visible here. Checks run in
# the order listed in check_name.
valid_checks <- c(
  "check_convergence",
  "aic",
  "mohns_rho",
  "compare_wham_models",
  "plot_wham_output",
  "self_test"
)
unknown_checks <- setdiff(check_name, valid_checks)
if (length(unknown_checks)) {
  stop("Unknown check_name: ", paste(unknown_checks, collapse = ", "))
}

if ("check_convergence" %in% check_name) {
  wham::check_convergence(fit, ret = TRUE)
}

if ("aic" %in% check_name) {
  wham::aic(fit)
}

if ("mohns_rho" %in% check_name) {
  if (is.null(fit$peels)) stop("This fit has no retrospective peels.")
  wham::mohns_rho(fit)
}

if ("compare_wham_models" %in% check_name) {
  wham::compare_wham_models(
    list(model = fit),
    do.table = TRUE,
    do.plot = FALSE,
    fdir = output_dir,
    table.opts = list(
      calc.rho = FALSE,
      calc.aic = TRUE,
      print = FALSE,
      save.csv = FALSE
    )
  )
}

if ("plot_wham_output" %in% check_name) {
  wham::plot_wham_output(
    fit,
    dir.main = file.path(output_dir, "plot_wham_output"),
    # PNG/PDF output also creates LaTeX parameter tables. HTML avoids the
    # LaTeX dependency while still generating the PNG plots.
    out.type = "html",
    plot.opts = list(browse = FALSE)
  )
}

if ("self_test" %in% check_name) {
  fit_rds <- file.path(output_dir, "fit.RDS")
  saveRDS(fit, fit_rds)
  wham::self_test(
    fit_RDS = fit_rds,
    n = 1,
    seeds = 8675,
    do_parallel = FALSE
  )
}

# Run this only after the selected checks have been run successfully.
if (run_refit) {
  # The current refit_wham default uses the stored input for compatibility.
  refit <- wham::refit_wham(fit)
  refit
}
