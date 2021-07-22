#' Return the posterior summary R^2 diagnostic 
#' 
#' Returns posterior samples of the R^2 diagnostic for psoterior summaries
#'
#' @param possum_fit A posterior summary object (e.g. returned by additive_summary)
#'
#' @return A vector of posterior samples
#' @export
#'
#' @examples
get_additive_summary_Rsq = function(possum_fit) {
  possum_fit$summaryRsq
}

#' Return posterior summaries of categorical variables
#' 
#' Return posterior samples of the coefficients on factor variables in a generalized additive posterior summary 
#'
#' @param possum_fit A posterior summary object (e.g. as returned by additive_summary)
#'
#' @return A matrix of posterior samples, with one column per coefficient
#' @export
#'
#' @examples
get_additive_factor_posterior = function(possum_fit) {
  possum_fit$gamFactorDf
}