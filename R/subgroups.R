##' Aggregate posterior draws into subgroups by a factor
##'
##' Aggregate posterior draws into subgroups by a factor and return a new
##' matrix of posterior draws for the subgroup (weighted) average
##' @title subgroup_average_posterior
##' @param posterior_samples A matrix of posterior samples, with MCMC draws down the rows
##' @param groups A factor defining the desired subgroups. Must have length ncol(posterior_samples)
##' @param weights Optional weights, for computing a weighted subgroup average
##' @return A matrix of size (number of MCMC samples) x (number of subgroups) 
##' containing posterior draws of the subgroup averages
##' 
subgroup_average_posterior = function(posterior_samples, groups, weights = NULL) {
  
  gpind = fastDummies::dummy_cols(data.frame(group = factor(groups)),
                                  remove_selected_columns = TRUE)
  gpind = as.matrix(gpind)
  
  colnames(gpind) = make.names(substring(colnames(gpind), first=7))
  
  if(!is.null(weights)) gpind = gpind*weights
  out = t(t(posterior_samples %*% gpind)/colSums(gpind))
  colnames(out) = colnames(gpind)
  
  return(data.frame(out))
  
}

##' Aggregate posterior draws by averaging
##'
##' Aggregate posterior draws by averaging, re (weighted) average
##' @title subgroup_average_posterior
##'
##' @param posterior_samples A matrix of posterior samples, with MCMC draws down the rows
##' @param weights Optional weights, for computing a weighted subgroup average
##'
##' @return
average_posterior = function(posterior_samples, weights = NULL) {
  
  groups = rep("ATE", ncol(posterior_samples))
  
  subgroup_average_posterior(posterior_samples, groups)
  
}

#' Difference_posterior
#'
#' @param posterior_samples Set of posterior samples
#' @param ref_col Differences are taken as (each column) - (ref_col)
#' @param remove_ref Remove the column of zeros corresponding to the reference column
#'
#' @return
#' @export
#'
#' @examples
difference_posterior = function(posterior_samples, ref_col = 1, remove_ref = TRUE) {
  diff_samples = posterior_samples - posterior_samples[,ref_col]
  if(remove_ref) diff_samples = diff_samples[,-ref_col]
  diff_samples
}