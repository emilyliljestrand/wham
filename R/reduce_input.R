#' Reduce the years of the model
#'
#' Internal function called by \code{\link{fit_peel}} for \emph{i} in 1--\code{n.peels}.
#' Creates the input for the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
#' Reduces the year dimension of all data, parameters, maps, so that one can project correctly from the peeled model.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}}). NOT from prepare_projection or project_wham
#' @param years_peeled which of input$years to peel from the model input
#' @param retro (T/F) whether this is for a retro peel (Default = TRUE)
#'
reduce_input <- function(input,years_peeled, retro = TRUE){
	
	if(input$data$n_years_proj>0) stop("cannot run wham:::reduce_input on an input for a projection")
		if (!length(years_peeled) || anyNA(years_peeled) ||
			anyDuplicated(years_peeled) || !all(years_peeled %in% input$years)) {
			stop("years_peeled must identify distinct model years.")
		}
		if (!setequal(years_peeled, tail(input$years, length(years_peeled))) ||
			length(input$years) - length(years_peeled) < 2L) {
			stop("reduce_input requires terminal-year peels retaining at least two model years.")
		}
	new <- input
	data <- new$data
	par <- new$par
	map <- new$map
	ind <- which(!(input$years %in% years_peeled))
	years <- input$years[ind]
	n_years <- length(years)
	new$years <- new$years_full <- years


	data$n_years_model <- n_years
	data$years_use <- data$years_use[ind]
	data$XSPR_R_avg_yrs <- data$XSPR_R_avg_yrs[which(data$XSPR_R_avg_yrs %in% (ind-1))]
	if(length(data$XSPR_R_avg_yrs)<1) {
		warning("Years specified to use for average recruitment in SPR-based reference points are not included in years of reduced input so all years will be used.")
		data$XSPR_R_avg_yrs <- 1:n_years - 1
	}
	# data$avg_years_ind <- data$n_years_model - input$data$n_years_model + data$avg_years_ind#shift which years to average back given the number of years to peel.
	# if(all(data$avg_years_ind<0)) data$avg_years_ind <- 1:data$n_years_model - 1 #shouldn't really ever happen
	# data$avg_years_ind <- data$avg_years_ind[which(data$avg_years_ind>=0)] #reduce if number of years used is more than the number available in the peel.

	avg_years_ind_L <- avg_years_ind_static_L <- matrix(0, data$n_years_model+1, data$n_regions) 
	for(r in 1:data$n_regions){

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_L[2:(input$data$n_years_model+1),r] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_L[1:(sum(temp>=0)+1),r] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_L[2:(input$data$n_years_model+1),r] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_L[1:(sum(temp>=0)+1),r] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

	}
	data$avg_years_ind_L <- avg_years_ind_L
	data$avg_years_ind_static_L <- avg_years_ind_static_L

	avg_years_ind_sel <- avg_years_ind_waacatch <- avg_years_ind_static_sel <- avg_years_ind_static_waacatch <- matrix(0, data$n_years_model+1, data$n_fleets) 
	for(f in 1:data$n_fleets){
		
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_sel[2:(input$data$n_years_model+1),f] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_sel[1:(sum(temp>=0)+1),f] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
		
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_waacatch[2:(input$data$n_years_model+1),f] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_waacatch[1:(sum(temp>=0)+1),f] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_sel[2:(input$data$n_years_model+1),f] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_sel[1:(sum(temp>=0)+1),f] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
		
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_waacatch[2:(input$data$n_years_model+1),f] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_waacatch[1:(sum(temp>=0)+1),f] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
	}
	data$avg_years_ind_sel <- avg_years_ind_sel
	data$avg_years_ind_waacatch <- avg_years_ind_waacatch
	data$avg_years_ind_static_sel <- avg_years_ind_static_sel
	data$avg_years_ind_static_waacatch <- avg_years_ind_static_waacatch

	avg_years_ind_waassb <- avg_years_ind_mat <- avg_years_ind_move <- avg_years_ind_M <- array(0, dim = c(data$n_stocks, data$n_regions, data$n_years_model+1)) 
	avg_years_ind_static_waassb <- avg_years_ind_static_mat <- avg_years_ind_static_move <- avg_years_ind_static_M <- avg_years_ind_static_SRR <- avg_years_ind_waassb
	avg_years_ind_NAA <- array(0, dim = c(data$n_stocks, data$n_regions, data$n_ages, data$n_years_model+1)) 
	for(s in 1:data$n_stocks) for(r in 1:data$n_regions){
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_waassb[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_waassb[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
		
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_mat[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_mat[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_move[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_move[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_M[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_M[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		for(a in 1:data$n_ages){
			temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_NAA[s,r,a,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
			if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
			avg_years_ind_NAA[s,r,a,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
		}

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_waassb[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_waassb[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
		
		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_mat[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_mat[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_move[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_move[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_M[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_M[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.

		temp <- data$n_years_model - input$data$n_years_model + data$avg_years_ind_static_SRR[s,r,2:(input$data$n_years_model+1)] #shift which years to average back given the number of years to peel.
		if(all(temp<0)) temp <- 1:data$n_years_model - 1 #shouldn't really ever happen
		avg_years_ind_static_SRR[s,r,1:(sum(temp>=0)+1)] <- c(sum(temp>=0),temp[which(temp>=0)]) #reduce if number of years used is more than the number available in the peel.
	}
	data$avg_years_ind_waassb <- avg_years_ind_waassb
	data$avg_years_ind_mat <- avg_years_ind_mat
	data$avg_years_ind_M <- avg_years_ind_M
	data$avg_years_ind_move <- avg_years_ind_move
	data$avg_years_ind_NAA <- avg_years_ind_NAA

	data$avg_years_ind_static_waassb <- avg_years_ind_static_waassb
	data$avg_years_ind_static_mat <- avg_years_ind_static_mat
	data$avg_years_ind_static_M <- avg_years_ind_static_M
	data$avg_years_ind_static_move <- avg_years_ind_static_move
	data$avg_years_ind_static_SRR <- avg_years_ind_static_SRR


	# data$avg_years_ind_static <- data$avg_years_ind

	data$which_F_age <- data$which_F_age[ind]
	data$FXSPR_init <- data$FXSPR_init[ind]
	data$FMSY_init <- data$FMSY_init[ind]

	data$fracyr_SSB <- data$fracyr_SSB[ind,,drop=F]
	data$mature <- data$mature[,ind,,drop=F]
	data$catch_Neff <- data$catch_Neff[ind,,drop=F]
	data$agg_catch <- data$agg_catch[ind,,drop=F]
	data$agg_catch_sigma <- data$agg_catch_sigma[ind,,drop=F]
	data$catch_paa <- data$catch_paa[,ind,,drop=F]
	data$use_agg_catch <- data$use_agg_catch[ind,,drop=F]
	data$use_catch_paa <- data$use_catch_paa[ind,,drop=F]
	data$selblock_pointer_fleets <- data$selblock_pointer_fleets[ind,,drop=F]
	data$index_Neff <- data$index_Neff[ind,,drop=F]
	data$agg_index_sigma <- data$agg_index_sigma[ind,,drop=F]
	data$agg_indices <- data$agg_indices[ind,,drop=F]
	data$index_paa <- data$index_paa[,ind,,drop=F]
	data$use_indices <- data$use_indices[ind,,drop=F]
	data$use_index_paa <- data$use_index_paa[ind,,drop=F]
	data$selblock_pointer_indices <- data$selblock_pointer_indices[ind,,drop=F]
	data$fracyr_indices <- data$fracyr_indices[ind,,drop=F]
	data$waa <- data$waa[,ind,,drop=F]
	data$selblock_years <- data$selblock_years[ind,,drop=F]
	data$n_years_selblocks <- apply(data$selblock_years,2,sum)

	data$keep_C <- data$keep_C[ind,,drop=F]
	data$keep_I <- data$keep_I[ind,,drop=F]
	data$keep_Cpaa <- data$keep_Cpaa[,ind,,drop=F]
	data$keep_Ipaa <- data$keep_Ipaa[,ind,,drop=F]

	if(!is.null(map$q_re)) {
		tmp <- par$q_re
		tmp[] <- as.integer(map$q_re)
		map$q_re <- factor(tmp[ind,])	
	}
	par$q_re <- par$q_re[ind,,drop=F]

	if(!is.null(map$mu_re)) {
		tmp <- par$mu_re
		tmp[] <- as.integer(map$mu_re)
		map$mu_re <- factor(tmp[,,,ind,,])	
	}
	par$mu_re <- par$mu_re[,,,ind,,,drop=F]

	if(!is.null(map$M_re)) {
		tmp <- par$M_re
		tmp[] <- as.integer(map$M_re)
		map$M_re <- factor(tmp[,,ind,])	
	}
	par$M_re <- par$M_re[,,ind,,drop=F]
	
	if(!is.null(map$log_NAA)) {
		tmp <- par$log_NAA
		tmp[] <- as.integer(map$log_NAA)
		map$log_NAA <- factor(tmp[,,head(ind,n_years-1),])	
	}
	par$log_NAA<- par$log_NAA[,,head(ind,n_years-1),,drop = F] 
	
	if(!is.null(map$L_re)) {
		tmp <- par$L_re
		tmp[] <- as.integer(map$L_re)
		map$L_re <- factor(tmp[ind,])	
	}
	par$L_re <- par$L_re[ind,,drop=F]
	
	if(!is.null(map$selpars_re)){
		tmp <- par$selpars_re
		tmp[] <- as.integer(map$selpars_re)
		map$selpars_re <- factor(tmp[,ind,]) #n_selblocks x n_years x n_ages
  }
  par$selpars_re <- par$selpars_re[,ind,,drop=F]

		if (!is.null(map$F_pars)) {
			tmp <- matrix(as.character(map$F_pars), nrow = NROW(par$F_pars))
			map$F_pars <- factor(tmp[ind, , drop = FALSE])
		}
		par$F_pars <- par$F_pars[ind, , drop = FALSE]

		# A terminal selectivity block can disappear completely in a peel.
		empty_blocks <- which(data$n_years_selblocks == 0L)
		if (length(empty_blocks)) {
			data$selblock_models_re[empty_blocks] <- 0L
			for (field in c("logit_selpars", "sel_repars", "selpars_re")) {
				mapping <- map[[field]]
				if (is.null(mapping)) mapping <- factor(seq_along(par[[field]]))
				mapping <- array(as.character(mapping), dim(par[[field]]))
				if (length(dim(mapping)) == 2L) mapping[empty_blocks, ] <- NA
				else mapping[empty_blocks, , ] <- NA
				map[[field]] <- factor(mapping)
			}
		}

  dims <- dim(par$q_re)
  dims[2] <- data$n_years_model + data$n_years_proj
  tmp <- matrix(0, nrow = data$n_years_model + data$n_years_proj, ncol = NCOL(par$q_re))
  tmp[1:data$n_years_model,] <- par$q_re[1:data$n_years_model,]
  par$q_re <- tmp

	new$par <- par
	new$map <- map
	new$data <- data

		new <- reduce_ecov_input(new, input)
		new$random <- Filter(function(field) {
			is.null(new$map[[field]]) || any(!is.na(new$map[[field]]))
		}, new$random)

	#reset the obs_vec
	new$by_pwi <- TRUE #needed to check age ranges of age comp
	new <- set_osa_obs(new)
	new$by_pwi <- NULL
	return(new)

}

reduce_ecov_input <- function(new, input) {
  if (!any(input$data$Ecov_model > 0) && is.null(input$options$ecov)) return(new)
  years <- input$years_Ecov
  if (!is.numeric(years) || length(years) != input$data$n_years_Ecov ||
      any(!is.finite(years)) || any(diff(years) != 1)) {
    stop("Cannot safely peel environmental covariates: years_Ecov must match the consecutive Ecov data years.")
  }
  keep <- which(years <= tail(new$years, 1L))
  if (!length(keep)) stop("No environmental covariate years remain in this peel.")
  new$years_Ecov <- years[keep]
  new$data$n_years_Ecov <- length(keep)
  new$data$years_use_Ecov <- input$data$years_use_Ecov[input$data$years_use_Ecov %in% (keep - 1L)]
  for (field in c("Ecov_obs", "Ecov_use_obs")) {
    new$data[[field]] <- input$data[[field]][keep, , drop = FALSE]
  }
  for (field in c("Ecov_re", "Ecov_obs_logsigma", "Ecov_obs_logsigma_re")) {
    # Old inputs can have unused placeholder sigma RE matrices.
    if (NROW(input$par[[field]]) != length(years)) {
      if (field == "Ecov_obs_logsigma_re" && !any(input$data$Ecov_obs_sigma_opt == 4L)) next
      stop("Cannot safely peel environmental covariates: invalid ", field, " dimensions.")
    }
    new$par[[field]] <- input$par[[field]][keep, , drop = FALSE]
    if (!is.null(input$map[[field]])) {
      mapping <- matrix(as.character(input$map[[field]]), nrow = length(years))
      new$map[[field]] <- factor(mapping[keep, , drop = FALSE])
    }
  }
  for (effect in c("R", "M", "q", "mu")) {
    field <- paste0("ind_Ecov_out_end_", effect)
    new$data[[field]] <- input$data[[field]] - (length(input$years) - length(new$years))
  }
  new$data$avg_years_Ecov <- input$data$avg_years_Ecov[input$data$avg_years_Ecov %in% (keep - 1L)]
  ecov <- input$options$ecov
  if (!is.null(ecov)) {
    rows <- which(ecov$year <= tail(new$years, 1L))
    old_n <- length(ecov$year)
    ecov$year <- ecov$year[rows]
    ecov$mean <- as.matrix(ecov$mean)[rows, , drop = FALSE]
    ecov$use_obs <- as.matrix(ecov$use_obs)[rows, , drop = FALSE]
    subset_sigma <- function(x) {
      if (is.matrix(x) && nrow(x) == old_n) return(x[rows, , drop = FALSE])
      if (is.numeric(x) && input$data$n_Ecov == 1L && length(x) == old_n) return(x[rows])
      x
    }
    if (is.list(ecov$logsigma)) ecov$logsigma[[1L]] <- subset_sigma(ecov$logsigma[[1L]])
    else ecov$logsigma <- subset_sigma(ecov$logsigma)
    new$options$ecov <- ecov
  }
  new
}