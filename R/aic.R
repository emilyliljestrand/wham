#' Calculate Akaike's Information Criterion for a WHAM model
#'
#' Calculates AIC from a fitted WHAM model's minimized negative log-likelihood
#' and number of estimated parameters.
#'
#' @param mod A fitted WHAM model object returned by [fit_wham()].
#'
#' @return A numeric AIC value.
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- fit_wham(input)
#' aic(mod)
#' }
aic <- function(mod) {
  2 * (mod$opt$obj + length(mod$opt$par))
}