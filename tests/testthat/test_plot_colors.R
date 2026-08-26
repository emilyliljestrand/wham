test_that("WHAM plot color modes return the expected palettes", {
  old.options <- options(wham.colors = "default")
  on.exit(options(old.options), add = TRUE)

  expect_equal(wham:::wham_palette(4), grDevices::hcl.colors(4, palette = "Blues 2", rev = FALSE))

  options(wham.colors = "classic")
  expect_equal(wham:::wham_palette(4), grDevices::hcl.colors(4, palette = "Blues 2", rev = FALSE))

  options(wham.colors = "bw")
  bw <- grDevices::col2rgb(wham:::wham_palette(4))
  expect_equal(bw[1, ], bw[2, ])
  expect_equal(bw[2, ], bw[3, ])

  options(wham.colors = "muted")
  expect_length(wham:::wham_palette(4), 4)
  expect_false(identical(wham:::wham_palette(4), grDevices::hcl.colors(4, palette = "Blues 2", rev = FALSE)))

  options(wham.colors = "viridis")
  expect_equal(wham:::wham_palette(4), viridisLite::viridis(4))
})