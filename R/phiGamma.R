
##' Inflation of sigma
##'
##' .. content for \details{} ..
##' @title Inflation of stdev
##' @param y 
##' @param gamma 
##' @param sigma 
##' @return 
##' @author Spencer Woody
phiGamma <- function(y, gamma, sigma) {

  if (is.matrix(gamma)) {
    y <- matrix(rep(y, ncol(gamma)), ncol = ncol(gamma))

    sqrt(colMeans((y - gamma)^2)) / sigma - 1
  } else {
    sqrt(mean((y - gamma)^2)) / sigma - 1
  }
  

}
