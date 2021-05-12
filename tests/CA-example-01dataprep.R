
library(here)
library(ggplot2)
library(stringr)
library(dplyr)
library(rpart)
library(rpart.plot)
## library(ggmap)
library(latex2exp)
library(splines)
## library(laGP)

set.seed(1)

alpha <- 0.05

###############################################################################
                                        #             Read in data            #
###############################################################################

housing <- read.csv(here("data/calif_penn_2011.csv")) %>%
  na.omit()

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
