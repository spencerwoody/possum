##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Rsquared
##' @return
##' @author Spencer Woody
##'
##' @export
rsq <- function(y, yhat) {
  ssr <- sum((y - yhat)^2)
  sst <- sum((y - mean(y))^2)
  1 - ssr / sst
}
