
library(horseshoe)
library(ggplot2)
library(cowplot)
library(dplyr)
library(lars)
library(stringr)
library(tidyr)
library(latex2exp)
library(dbarts)
library(mgcv)

theme_set(theme_minimal_grid())

devtools::load_all()

###############################################################################
                                        #        Load in the cali data        #
###############################################################################

data(calif)

###############################################################################
                                        #          Train / test split         #
###############################################################################

## Add county column
calif <- housing %>%
  filter(STATEFP == 6) %>%
  mutate(County = str_extract(GEO.display.label, regex(", .* County")) %>%
           str_sub(3, -1))

calif %>% glimpse()

###############################################################################
                                        #       Response and covariates       #
###############################################################################

y <- calif %>%
  dplyr::pull(Median_house_value) %>%
  log()

califDf <- calif %>%
  dplyr::select(Median_household_income,
                POPULATION,
                Median_rooms,
                LONGITUDE,
                LATITUDE) %>%
  mutate(Median_household_income = log(Median_household_income),
         POPULATION = log(POPULATION))

terms <- califDf %>% names()

x <- califDf %>% as.matrix()

N <- nrow(x)

d <- ncol(x)

## source("tests/CA-example-01dataprep.R")

colnames(x)

###############################################################################
                                        #               run bart              #
###############################################################################

rsq <- function(y, yhat) {
  ssr <- sum((y - yhat)^2)
  sst <- sum((y - mean(y))^2)
  1 - ssr / sst
}

mybart <- bart(x, y, ndpost=1000, nskip=1000, nchain = 1)

## mybart <- bart2(x, y, n.samples=1000, n.burn=10000, n.chains=4)

## yhat <- mybart$yhat.train.mean
## yhatSamples <- mybart$yhat.train %>% apply(3, rbind) %>% t()

names(mybart)

yhat <- mybart$yhat.train.mean
yhatSamples <- t(mybart$yhat.train)

rsq(y, yhat)

plot(y, yhat)
plot(y, y-yhat)

dim(yhatSamples)

###############################################################################
                                        #               Use the               #
###############################################################################



myfactor <- sample(c("apples", "oranges"), nrow(califDf), replace=TRUE)
myfactor2 <- sample(c("beef", "chicken"), nrow(califDf), replace=TRUE)

califDf2 <- califDf %>%
  mutate(fruit = myfactor,
         fruit=factor(fruit)) %>%
  mutate(meat = myfactor2,
         meat=factor(meat)) %>%
  glimpse()

myform <- yhat ~
  s(Median_household_income) +
  s(POPULATION) +
  s(Median_rooms) +
  s(LONGITUDE) +
  s(LATITUDE)

myform <- yhat ~
  s(Median_household_income) +
  s(POPULATION) +
  s(Median_rooms) +
  s(LONGITUDE) +
  s(LATITUDE) +
  factor(fruit) +
  factor(meat)

myform <- yhat ~
  s(Median_household_income) +
  s(POPULATION) +
  s(Median_rooms) +
  (LONGITUDE) +
  s(LATITUDE) +
  factor(fruit) +
  factor(meat)

myform_lin <- yhat ~
                s(Median_household_income) +
                s(POPULATION) +
                s(Median_rooms) +
                s(LONGITUDE) +
                LATITUDE

args(gamProjectionFun)

mytry_tri <- gamProjectionFun(myform, yhatSamples, yhat, califDf2, verbose = TRUE)

names(mytry_tri)

mytry_tri$gamPlot

mytry_tri$gamTerm

mytry_tri$gamFactorDf %>% glimpse()

mytry_tri$gamFactorDf %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_violin(aes(level, post), alpha=0.5, fill = "grey90") +
  facet_wrap(~term, scales="free_x")

mytry_tri$gamFactorDf %>%
  ggplot()

mytry_tri <- gamProjectionFun(myform, yhatSamples, yhat, califDf, verbose = TRUE)

mytry_tri %>% names()












gamFit <- gam(myform,
               data = califDf2)


Xgam <- model.matrix(gamFit)
V <- vcov(gamFit, dispersion = 1)
Q <- crossprod(V, crossprod(Xgam, fhatmat))













mytry <- additiveSummary1(myform, yhatSamples, yhat, califDf)

mytry$gamPlot

mytry_lin <- additiveSummary1(myform_lin, yhatSamples, yhat, califDf)

mytry_lin$gamPlot

foobar <- califDf3 %>% select(terms, ran_factor) %>% glimpse() 

foobar

mytry_lin$gamFit %>% summary()
mytry_lin$gamFit %>% confint()

myfoo <- mytry_lin$gamFit %>% predict(type="terms")

plot(califDf$LATITUDE, myfoo[, 1])
plot(califDf$LONGITUDE, myfoo[, 5])

mytry_lin$gamFit %>% model.matrix() %>% colnames()

mytry_lin$gammaTerms[["LATITUDE"]] %>% dim()

foo <- mytry_lin$gammaTerms[["POPULATION"]]
bar <- califDf$POPULATION

foo <- mytry_lin$gammaTerms[["LONGITUDE"]]
bar <- califDf$LONGITUDE

100 * 99 / 2

50 * 49 / 2

myquants <- quantile(bar, probs = seq(0.025, 0.975, length.out = 100), type=3)

myinds <- which(bar %in% myquants)
myinds_len <- length(myquants)

###############################################################################
                                        #             Triangles...            #
###############################################################################

ivec <- rep(NA, myinds_len * (myinds_len - 1) / 2)
jvec <- rep(NA, myinds_len * (myinds_len - 1) / 2)
x_lo <- rep(NA, myinds_len * (myinds_len - 1) / 2)
x_hi <- rep(NA, myinds_len * (myinds_len - 1) / 2)
prob <- rep(NA, myinds_len * (myinds_len - 1) / 2)

idx <- 1

for (i in 1:(myinds_len - 1)) {
  for (j in (i+1):myinds_len) {
    ## mean(foo[myinds[i], ] < foo[myinds[j], ])
    ivec[idx] <- i
    jvec[idx] <- j
    x_lo[idx] <- bar[myinds[i]]
    x_hi[idx] <- bar[myinds[j]]
    prob[idx] <- mean(foo[myinds[j], ] > foo[myinds[i], ])
    idx <- idx+1
  }
}

triangle <- data.frame(i=ivec, j=jvec, prob=prob)

triangle %>%
  ggplot() +
  geom_raster(aes(i, j, fill = prob)) +
  scale_fill_gradient2(midpoint=0.5) +
  labs(x = "lower quantile", y = "upper quantile") +
  coord_equal()

bar[bar %in% myquants] %>% length()

quantile(1:7, type=7)
quantile(1:7, type=3)

bar[bar %in% myquants]

cut(bar, quantile(bar))

califDf_quant <- califDf %>%
  mutate(quantile = ntile(LATITUDE, 100))

names(mytry_lin$gammaTerms) <- terms

windsorize <- function(gamDf, windsor = 0.05) {}




###############################################################################
                                        #             GAM summary             #
###############################################################################

myform <- yhat ~
                s(Median_household_income) +
                s(POPULATION) +
                s(Median_rooms) +
                s(LONGITUDE) +
                (LATITUDE)

califDf3 <- califDf %>%
  mutate(ran = 2 * rbinom(n(), 1, 0.5) + rbinom(n(), 1, 0.5)) %>%
  mutate(ran = ran) %>% 
  mutate(ran_fact = factor(ran)) %>%
  glimpse()

is.factor(califDf3$ran)
is.factor(califDf3$ran_fact)

unique(califDf3$ran_fact)


gamFit3 <- gam(yhat ~
                 s(Median_household_income) +
                 s(POPULATION) +
                 s(Median_rooms) +
                 s(LONGITUDE) +
                 s(LATITUDE) +
                 ran_fact,
               data = califDf3)

mymat <- gamFit3 %>% model.matrix()

gamFit3 %>% summary()

colnames(mymat)[1:10]

glimpse(mymat)

mymat[1:10, 2:4]
califDf3$ran_fact[1:10]

mymat[, "ran_fact1"]

gamFit3 %>% model.matrix() %>% apply(2, function(x) length(unique(x)))

gamFit <- gam(myform,
              data = califDf)

colnames(model.matrix(gamFit))

## Model matrix and covariance matrix
Xs <- model.matrix(gamFit)
V <- vcov(gamFit, dispersion = 1)

## Point summary and posterior of projection
q <- crossprod(V, crossprod(Xs, yhat))
Q <- crossprod(V, crossprod(Xs, yhatSamples))

## Posterior of fitted values from projections
gammaSamples <- Xs %*% Q

gamma <- predict(gamFit, newdata = califDf, type = "response")
gammaTerms <- predict(gamFit, newdata = califDf, type = "terms",
                      se.fit = TRUE)

gammaTerms$fit

gamres <- as.numeric(yhat - as.numeric(gamma))

gammaTerms[[1]] %>% dim()
gammaTerms[[2]] %>% dim()

length(gammaTerms)

colnames(Xs)

gsub("\\([^()]*\\)", "", colnames(Xs))


gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", colnames(Xs2), perl=T)

