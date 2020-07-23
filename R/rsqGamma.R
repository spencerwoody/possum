
##' Description
##'
##' Details
##' @title A function for computing the summary R-squared
##' @param gamma fitted values of summary (vector or matrix)
##' @param fhat fitted values from regression function (vector or
##'   matrix)
##' @return
##' @author Spencer Woody
##' @export
rsqGamma <- function(gamma, fhat) {

  if (!is.matrix(gamma)) {
    gamma <- matrix(gamma, ncol = 1)
  }

  if (!is.matrix(fhat)) {
    fhat <- matrix(fhat, ncol = 1)
  }

  fhatBar <- colMeans(fhat)

  SST <- colSums(sweep(fhat, 2, fhatBar)^2)

  if (ncol(gamma) == 1) {
    gammaMat <- matrix(rep(as.numeric(gamma), ncol(fhat)),
                       ncol = ncol(fhat))

    SSR <- colSums((fhat - gammaMat)^2)
  } else {
    SSR <- colSums((fhat - gamma)^2)
  }

  1 - SSR / SST

}
