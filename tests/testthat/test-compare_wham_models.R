test_that("legacy fits can be compared with current refits", {
  for (model in c("Butterfishfit", "GBK_Cod_2024MTfit")) {
    original <- test_path("..", "..", "oldfitobjects", paste0(model, ".rds"))
    refit <- test_path(
      "..",
      "..",
      "oldfitobjects",
      "compare_fit_refit",
      model,
      paste0(model, "_refit.rds")
    )
    skip_if_not(
      all(file.exists(c(original, refit))),
      "Local original/refit pair is unavailable"
    )
    output <- tempfile(
      "compare-wham-",
      tmpdir = Sys.getenv("PA_SCRATCH_DIR", tempdir())
    )
    dir.create(output)
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    mods <- list(fit = readRDS(original), refit = readRDS(refit))
    result <- suppressMessages(compare_wham_models(
      mods,
      fdir = output,
      table.opts = list(calc.rho = FALSE, print = FALSE, save.csv = TRUE),
      plot.opts = list(browse = FALSE)
    ))
    # The Kobe plot is drawn directly rather than returned as a ggplot.
    expect_equal(
      vapply(result$g[1:9], inherits, logical(1), "ggplot"),
      rep(TRUE, 9)
    )
    expect_length(
      list.files(file.path(output, "compare_png"), pattern = "\\.png$"),
      10L
    )
    expect_equal(unname(result$daic), c(0, 0))

    mods$refit$env$data$do_SPR_BRPs <- 0L
    result <- suppressMessages(compare_wham_models(
      mods,
      do.table = FALSE,
      fdir = output,
      plot.opts = list(which = 8:10, browse = FALSE)
    ))
    expect_equal(sum(vapply(result$g, is.null, logical(1))), 10L)
  }
})

test_that("fits without sdreport can be compared with saved refits", {
  for (model in c("tfit_d_d", "tfit_dm_dm_vast_constant", "tfit_dm_dm")) {
    original <- test_path("..", "..", "oldfitobjects", paste0(model, ".rds"))
    refit <- test_path(
      "..",
      "..",
      "oldfitobjects",
      "compare_fit_refit",
      model,
      paste0(model, "_refit.rds")
    )
    skip_if_not(
      all(file.exists(c(original, refit))),
      "Local original/refit pair is unavailable"
    )
    output <- tempfile(
      "compare-wham-",
      tmpdir = Sys.getenv("PA_SCRATCH_DIR", tempdir())
    )
    dir.create(output)
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    mods <- list(fit = readRDS(original), refit = readRDS(refit))
    warnings <- character()
    result <- withCallingHandlers(
      suppressMessages(compare_wham_models(
        mods,
        fdir = output,
        table.opts = list(calc.rho = FALSE, print = FALSE, save.csv = TRUE),
        plot.opts = list(browse = FALSE)
      )),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_equal(
      warnings,
      "Not plotting SPR-based reference points and status because at least one model does not have them reported."
    )
    expect_equal(
      vapply(result$g[1:7], inherits, logical(1), "ggplot"),
      rep(TRUE, 7)
    )
    expect_equal(vapply(result$g[8:10], is.null, logical(1)), rep(TRUE, 3))
    expect_length(list.files(file.path(output, "compare_png"), "\\.png$"), 7L)
    expect_equal(nrow(result$g[[1]]$data), 6L * length(mods$fit$years_full))
    expect_equal(
      result$g[[2]]$data$val,
      rep(0, 6L * length(mods$fit$years_full))
    )
    expect_equal(nrow(result$tab), 2L)
    expect_equal(length(result$aic), 2L)
  }
})
