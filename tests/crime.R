
library(horseshoe)
library(ggplot2)
library(cowplot)
library(dplyr)
library(lars)
library(stringr)
library(tidyr)

theme_set(theme_minimal_grid())

devtools::load_all()

###############################################################################
                                        #         Prepare UScrime data        #
###############################################################################

data(UScrime, package = "MASS")

y <- UScrime %>%
  pull(y) %>%
  scale()

## Covariates
X <- UScrime[, -which(colnames(UScrime) == "y")] %>%
  as.matrix()

## log-transform and scale the data
X[, -2] <- log(X[, -2])

X <- scale(X)

varnames <- colnames(X)

varnamesDf <- data.frame(
  Var1 = seq_along(varnames),
  varname = varnames
)

###############################################################################
                                        #       Fit model with horseshoe     #
###############################################################################

fit <- horseshoe(y, X, method.tau="halfCauchy",
                 method.sigma="Jeffreys", nmc=45000, thin=5)

betaSamples <- fit$BetaSamples

## Posterior samples of beta; each covariate is one row
betaSamples <- fit$BetaSamples
sigma2Samples <- fit$Sigma2Samples

devtools::load_all()

sparseProj <- sparseLinearProj(X, y,
                               betaSamples, sigma2Samples,
                               varnames = varnames)

names(sparseProj)

df <- sparseProj$summaryDf %>%
  group_by(modelSize) %>%
  summarize(phi_gamma_lo = quantile(phi_gamma, 0.05),
            phi_gamma_hi = quantile(phi_gamma, 0.05),
            phi_gamma_mid = mean(phi_gamma),
            rsq_gamma_lo = quantile(rsq_gamma, 0.05),
            rsq_gamma_hi = quantile(rsq_gamma, 0.05),
            rsq_gamma_mid = mean(rsq_gamma))

df %>%
  ggplot(aes(modelSize, rsq_gamma_mid)) +
  geom_point() 

sparseProj$summaryDf %>% head()
sparseProj$summaryDfLong %>% head()


###############################################################################
                                        #    Nonparametric linear summaries   #
###############################################################################

library(dbarts)

cbart <- bart(X, y)

plot(cbart$yhat.train.mean, y)

yhatmat <- t(cbart$yhat.train)
fhatmat <- t(cbart$yhat.train)

devtools::load_all()

mytry <- sparseLinearSummary(X, fhatmat, y, sigma2Samples)

names(mytry)

mytry$modelSize

df2 <- mytry$summaryDf %>%
  group_by(modelSize) %>%
  summarize(phi_gamma_lo = quantile(phi_gamma, 0.05),
            phi_gamma_hi = quantile(phi_gamma, 0.95),
            phi_gamma_mid = mean(phi_gamma),
            rsq_gamma_lo = quantile(rsq_gamma, 0.05),
            rsq_gamma_hi = quantile(rsq_gamma, 0.95),
            rsq_gamma_mid = mean(rsq_gamma))

df2 %>%
  ggplot(aes(modelSize, rsq_gamma_mid)) +
  geom_linerange(aes(ymin = rsq_gamma_lo,
                     ymax = rsq_gamma_hi)) + 
  geom_point() +
  labs()

head(m)

mytry$rsqGammaList
