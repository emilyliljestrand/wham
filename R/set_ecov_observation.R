# Configure observation sigma parameters and maps; return updated input, ecov, data, par, and map.
# Preserve the original logsigma specification and normalize process_model across covariates.
set_ecov_observation <- function(input, ecov, data, par, map, n_Ecov_obs) {
  #default NAs for parameter matrix of observation standard deviations
  par$Ecov_obs_logsigma <- matrix(NA, n_Ecov_obs, data$n_Ecov)
  map$Ecov_obs_logsigma <- matrix(NA, nrow = n_Ecov_obs, ncol = data$n_Ecov) # turn off estimation
  par$Ecov_obs_sigma_par <- matrix(0, nrow = 2, ncol = data$n_Ecov)
  map$Ecov_obs_sigma_par <- matrix(NA, nrow = 2, ncol = data$n_Ecov)

  logsigma_input <- ecov$logsigma # retain fixed values as well as estimation options for rebuilds
  logsig_more <- list()
  if (is.list(ecov$logsigma)) {
    n <- length(ecov$logsigma)
    if (n > 1) {
      logsig_more <- ecov$logsigma[[2]]
    } #further elements than the second are ignored.
    ecov$logsigma <- ecov$logsigma[[1]]
  }

  if (class(ecov$logsigma)[1] == "matrix") {
    par$Ecov_obs_logsigma <- ecov$logsigma
    if (identical(dim(par$Ecov_obs_logsigma), c(1L, data$n_Ecov))) {
      par$Ecov_obs_logsigma <- matrix(
        rep(par$Ecov_obs_logsigma, each = n_Ecov_obs),
        ncol = data$n_Ecov
      )
    }
    if (!identical(dim(par$Ecov_obs_logsigma), dim(data$Ecov_obs))) {
      stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
  }
  if (class(ecov$logsigma)[1] == 'numeric') {
    #data$Ecov_obs_sigma_opt[] <- 1 #defined above
    input$log$ecov <- c(
      input$log$ecov,
      "ecov$logsigma is numeric. Coercing to a matrix... \n"
    )
    if (length(ecov$logsigma) == data$n_Ecov) {
      par$Ecov_obs_logsigma <- matrix(
        rep(ecov$logsigma, each = n_Ecov_obs),
        ncol = data$n_Ecov
      )
    }
    if (length(ecov$logsigma) == n_Ecov_obs & data$n_Ecov == 1) {
      par$Ecov_obs_logsigma <- matrix(ecov$logsigma, ncol = 1)
    }
    if (
      length(ecov$logsigma) != data$n_Ecov & length(ecov$logsigma) != n_Ecov_obs
    ) {
      stop(
        "ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations"
      )
    }
  }

  #set up and check length of ecov$process_model
  if (length(ecov$process_model) == 1) {
    ecov$process_model <- rep(ecov$process_model, data$n_Ecov)
  } #use the single value for all Ecovs
  if (length(ecov$process_model) != data$n_Ecov) {
    stop("length of ecov$process_model must be either 1 or the number of Ecovs")
  }

  #now over write ecov$logsigma with logsig_more if available because the fixed obs var matrices have been defined
  if (length(logsig_more)) {
    ecov$logsigma <- logsig_more
  }

  if (class(ecov$logsigma)[1] == 'character') {
    #check that estimation options are right
    if (!all(ecov$logsigma %in% c("est_1", "est_re", NA))) {
      stop(
        "ecov$logsigma or ecov$logsigma[[2]] is character and must be NA (do not estimate), 'est_1' (single variance parameter), or 'est_re' (iid re annual variance parameters)"
      )
    }
    if (length(ecov$logsigma) == 1) {
      ecov$logsigma <- rep(ecov$logsigma, data$n_Ecov)
    } #use the single value for all Ecovs
    #check length of estimation options
    if (length(ecov$logsigma) != data$n_Ecov) {
      stop(
        "length of ecov$logsigma when character must be either 1 or the number of Ecovs"
      )
    }

    for (i in 1:data$n_Ecov) {
      if (!is.na(ecov$logsigma[i])) {
        if (ecov$logsigma[i] == 'est_1') {
          # estimate 1 Ecov obs sigma for each Ecov
          data$Ecov_obs_sigma_opt[i] <- 2
          par$Ecov_obs_logsigma[, i] <- -1.3
          map$Ecov_obs_logsigma[, i] <- i
        }
      }
      if (!is.na(ecov$logsigma[i])) {
        if (ecov$logsigma[i] == 'est_re') {
          data$Ecov_obs_sigma_opt[i] <- 4
          map$Ecov_obs_logsigma[, i] <- NA # turn off estimation of fixed effects
          par$Ecov_obs_sigma_par[, i] <- c(-1.3, -2.3) # random effect pars
          map$Ecov_obs_sigma_par[, i] <- max(
            0,
            map$Ecov_obs_sigma_par,
            na.rm = T
          ) +
            1:2
          par$Ecov_obs_logsigma[
            !is.finite(par$Ecov_obs_logsigma[, i]),
            i
          ] <- par$Ecov_obs_sigma_par[1, i]
        }
      }
    }
  }
  ecov$logsigma <- logsigma_input

  list(input = input, ecov = ecov, data = data, par = par, map = map)
}

# Allocate annual observation-sigma random effects using the final padded data dimensions.
# Return updated parameters and maps, leaving fixed-sigma columns inactive.
set_ecov_observation_re <- function(data, par, map) {
  # Allocate annual observation-sigma random effects only after padding is complete.
  par$Ecov_obs_logsigma_re <- matrix(0, data$n_years_Ecov, data$n_Ecov)
  map$Ecov_obs_logsigma_re <- matrix(NA, data$n_years_Ecov, data$n_Ecov)
  for (i in seq_len(data$n_Ecov)) {
    if (data$Ecov_obs_sigma_opt[i] == 4) {
      map$Ecov_obs_logsigma_re[, i] <- max(
        0,
        map$Ecov_obs_logsigma_re,
        na.rm = TRUE
      ) +
        seq_len(data$n_years_Ecov)
      par$Ecov_obs_logsigma_re[, i] <- par$Ecov_obs_logsigma[, i]
    }
  }

  list(par = par, map = map)
}

# Validate consecutive observation years and their count before any padding occurs.
validate_ecov_observation_years <- function(years, n_Ecov_obs) {
  if (
    !is.numeric(years) ||
      !is.null(dim(years)) ||
      !length(years) ||
      any(!is.finite(years)) ||
      any(diff(years) != 1)
  ) {
    stop(
      "ecov$year must be a nonempty, finite numeric vector of consecutive years"
    )
  }
  if (length(years) != n_Ecov_obs) {
    stop("ecov$year is not the same length as # rows in ecov$mean")
  }
  invisible(NULL)
}
