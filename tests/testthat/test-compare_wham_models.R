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
