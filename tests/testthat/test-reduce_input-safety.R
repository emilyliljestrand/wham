reduce_safety_input <- function(
  process = "ar1",
  sigma = "est_re",
  ecov_years = 2001:2003
) {
  ecov <- list(
    label = "temperature",
    year = ecov_years,
    mean = matrix(seq_along(ecov_years), ncol = 1),
    use_obs = matrix(1, length(ecov_years), 1),
    logsigma = sigma,
    process_model = process,
    recruitment_how = matrix("controlling-lag-1-linear", 1, 1),
    M_how = array("lag-0-linear", c(1, 1, 4, 1)),
    q_how = matrix("lag-2-linear", 1, 1)
  )
  suppressMessages(prepare_wham_input(
    NAA_re = list(sigma = "rec"),
    basic_info = list(years = 2000:2005, ages = 1:4),
    ecov = ecov
  ))
}

test_that("peels subset padded Ecov years for every process and sigma model", {
  for (process in c("rw", "ar1")) {
    for (sigma in list(-1.7, "est_1", "est_re")) {
      input <- reduce_safety_input(process, sigma)
      before <- input
      peeled <- reduce_input(input, 2004:2005)
      keep <- which(input$years_Ecov <= 2003)
      expect_identical(input, before)
      expect_equal(peeled$years, 2000:2003)
      expect_equal(peeled$years_full, peeled$years)
      expect_equal(peeled$years_Ecov, input$years_Ecov[keep])
      expect_equal(peeled$data$n_years_Ecov, length(keep))
      for (field in c("Ecov_obs", "Ecov_use_obs")) {
        expect_equal(
          peeled$data[[field]],
          input$data[[field]][keep, , drop = FALSE]
        )
      }
      for (field in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
        expect_equal(
          peeled$par[[field]],
          input$par[[field]][keep, , drop = FALSE]
        )
        expected <- matrix(
          as.character(input$map[[field]]),
          nrow = length(input$years_Ecov)
        )
        expect_equal(
          as.character(peeled$map[[field]]),
          as.vector(expected[keep, , drop = FALSE])
        )
      }
      expect_equal(
        peeled$data$years_use_Ecov,
        input$data$years_use_Ecov[input$data$years_use_Ecov %in% (keep - 1L)]
      )
      expect_equal(
        peeled$data$avg_years_Ecov,
        input$data$avg_years_Ecov[input$data$avg_years_Ecov %in% (keep - 1L)]
      )
      expect_identical(validate_input_safety(peeled), peeled)
    }
  }
})

test_that("Ecov MLEs, custom maps and random choices survive without stored options", {
  input <- reduce_safety_input()
  input$options$ecov <- NULL
  for (field in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
    input$par[[field]][] <- seq_along(input$par[[field]]) / 10 + 2
    input$map[[field]] <- factor(rep(
      c(NA, 7, 7, 9),
      length.out = length(input$par[[field]])
    ))
  }
  input$par$Ecov_process_pars[] <- c(0.4, -0.8, 0.2)
  input$par$Ecov_obs_sigma_par[] <- c(-0.6, 0.3)
  input$random <- c("Ecov_re", "Ecov_obs_logsigma_re", "F_pars")
  peeled <- reduce_input(input, 2004:2005)
  keep <- which(input$years_Ecov <= 2003)
  expect_null(peeled$options$ecov)
  expect_identical(peeled$random, input$random)
  for (field in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
    expect_equal(peeled$par[[field]], input$par[[field]][keep, , drop = FALSE])
    expect_equal(
      as.character(peeled$map[[field]]),
      as.character(input$map[[field]])[keep]
    )
  }
  for (field in c(
    "Ecov_process_pars",
    "Ecov_obs_sigma_par",
    "Ecov_beta_R",
    "Ecov_beta_M",
    "Ecov_beta_q"
  )) {
    expect_identical(peeled$par[[field]], input$par[[field]])
    expect_identical(peeled$map[[field]], input$map[[field]])
  }
  for (effect in c("R", "M", "q", "mu")) {
    start <- paste0("ind_Ecov_out_start_", effect)
    end <- paste0("ind_Ecov_out_end_", effect)
    expect_identical(peeled$data[[start]], input$data[[start]])
    expect_equal(peeled$data[[end]], input$data[[end]] - 2)
    if (any(input$data[[paste0("Ecov_how_", effect)]] > 0)) {
      expect_equal(
        peeled$data[[end]] - peeled$data[[start]] + 1,
        input$data[[end]] * 0 + 4
      )
    }
  }
  expect_identical(validate_input_safety(peeled), peeled)
})

test_that("stored numeric observation sigmas are subset without changing their form", {
  for (sigma in list(
    seq(-2, -1, length.out = 6),
    matrix(seq(-2, -1, length.out = 6), 6, 1),
    list(matrix(-1.4, 6, 1), "est_re"),
    -1.2,
    "est_1"
  )) {
    input <- reduce_safety_input(sigma = sigma, ecov_years = 2000:2005)
    peeled <- reduce_input(input, 2004:2005, retro = FALSE)
    expect_equal(peeled$options$ecov$year, 2000:2003)
    expected <- sigma
    if (is.list(sigma)) {
      expected[[1]] <- sigma[[1]][1:4, , drop = FALSE]
    } else if (is.matrix(sigma)) {
      expected <- sigma[1:4, , drop = FALSE]
    } else if (is.numeric(sigma) && length(sigma) == 6) {
      expected <- sigma[1:4]
    }
    expect_identical(peeled$options$ecov$logsigma, expected)
  }
})

test_that("terminal selectivity blocks lose all estimable maps and RE models", {
  input <- reduce_safety_input()
  input$data$selblock_pointer_fleets[, 1] <- c(1, 1, 1, 1, 3, 3)
  input <- set_selectivity(
    input,
    list(
      n_selblocks = 3,
      model = rep("logistic", 3),
      re = c("none", "none", "ar1_y")
    )
  )
  input$random <- c("selpars_re", "Ecov_re")
  peeled <- reduce_input(input, 2004:2005)
  expect_equal(peeled$data$n_years_selblocks[3], 0)
  expect_equal(peeled$data$selblock_models_re[3], 0)
  for (field in c("logit_selpars", "sel_repars", "selpars_re")) {
    mapping <- array(
      as.character(peeled$map[[field]]),
      dim(peeled$par[[field]])
    )
    removed <- if (length(dim(mapping)) == 2) mapping[3, ] else mapping[3, , ]
    expect_equal(all(is.na(removed)), TRUE)
    original <- array(as.character(input$map[[field]]), dim(input$par[[field]]))
    retained <- if (length(dim(mapping)) == 2) {
      mapping[1:2, ]
    } else {
      mapping[1:2, 1:4, ]
    }
    expected <- if (length(dim(original)) == 2) {
      original[1:2, ]
    } else {
      original[1:2, 1:4, ]
    }
    expect_equal(is.na(retained), is.na(expected))
    expect_equal(
      outer(as.vector(retained), as.vector(retained), "=="),
      outer(as.vector(expected), as.vector(expected), "==")
    )
  }
  expect_identical(peeled$random, "Ecov_re")
  expect_identical(validate_input_safety(peeled), peeled)
})

test_that("F maps preserve fixed and shared entries, including an absent map", {
  input <- reduce_safety_input()
  input$par$F_pars[] <- seq_along(input$par$F_pars) / 10
  input$map$F_pars <- factor(c(NA, 8, 8, 3, 9, 9))
  peeled <- reduce_input(input, 2004:2005)
  expect_equal(as.character(peeled$map$F_pars), c(NA, "8", "8", "3"))
  expect_equal(peeled$par$F_pars, input$par$F_pars[1:4, , drop = FALSE])
  input$map$F_pars <- NULL
  expect_null(reduce_input(input, 2004:2005)$map$F_pars)
})

test_that("invalid peel years and unsafe Ecov metadata fail informatively", {
  input <- reduce_safety_input()
  for (years in list(numeric(), NA, c(2005, 2005), 1999)) {
    error <- tryCatch(reduce_input(input, years), error = identity)
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "distinct model years", fixed = TRUE)
  }
  for (years in list(2002, c(2003, 2005), 2001:2005, 2000:2005)) {
    error <- tryCatch(reduce_input(input, years), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "terminal-year peels retaining at least two",
      fixed = TRUE
    )
  }
  projected <- input
  projected$data$n_years_proj <- 1
  error <- tryCatch(reduce_input(projected, 2005), error = identity)
  expect_match(conditionMessage(error), "projection", fixed = TRUE)
  for (years in list(
    NULL,
    c(1997, 1999),
    rep(2000, length(input$years_Ecov))
  )) {
    broken <- input
    broken$years_Ecov <- years
    error <- tryCatch(reduce_input(broken, 2005), error = identity)
    expect_s3_class(error, "error")
    expect_match(
      conditionMessage(error),
      "years_Ecov must match the consecutive",
      fixed = TRUE
    )
  }
  input$par$Ecov_re <- input$par$Ecov_re[-1, , drop = FALSE]
  error <- tryCatch(reduce_input(input, 2005), error = identity)
  expect_match(
    conditionMessage(error),
    "invalid Ecov_re dimensions",
    fixed = TRUE
  )
})
