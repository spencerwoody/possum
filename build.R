
## Compile attributes for Rcpp scripts
Rcpp::compileAttributes()
## RcppArmadillo::RcppArmadillo.package.skeleton("possum")
usethis::use_rcpp()
usethis::use_rcpp_armadillo()

##

## Create new R script


## DLL errors
pkgbuild::clean_dll()
pkgbuild::compile_dll()



## Create Rd files
devtools::document()

## Check that package can be build
devtools::check()

devtools::check_built()

## Install locally
devtools::install(upgrade = "always")

###############################################################################
                                        #                 Test                #
###############################################################################

## Load locally
devtools::load_all()

N <- 1000

x1 <- rnorm(N)
x2 <- rnorm(N)
y <- x1 + 2 * x2 + rnorm(N)

min(x1)
quantile(x1, 0, 1)
min(x2)

mygam <- gamWrap(y ~ s(x1) + s(x2))
