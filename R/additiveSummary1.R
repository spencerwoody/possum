##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title  Additive summary
##' @param gamcall 
##' @param yhatSamples 
##' @param yhat 
##' @param df 
##' @return 
##' @author Spencer Woody
additiveSummary1 <- function(gamcall,
                             fhatSamples,
                             fhat = rowMeans(fhatSamples),
                             df = NULL,
                             terms = NULL,
                             alpha = 0.05
                             ) {

  if (any(str_detect(gamcall, "factor\\("))) stop("Please make instead of using factor")

  factor_terms <- df %>%
    select_if(is.factor) %>%
    colnames()

  gamFit <- gam(gamcall, data = df)
  ## Model matrix and covariance matrix
  Xs <- model.matrix(gamFit)
  V <- vcov(gamFit, dispersion = 1)

  ## Point summary and posterior of projection
  q <- crossprod(V, crossprod(Xs, fhat))
  Q <- crossprod(V, crossprod(Xs, fhatSamples))

  ## Posterior of fitted values from projections
  gammaSamples <- Xs %*% Q

  gamma <- predict(gamFit, newdata = df, type = "response")
  gammaTerms <- predict(gamFit, newdata = df, type = "terms",
                        se.fit = TRUE)

  ## Terms
  terms_vec <- ifelse(colnames(Xs) %>% str_detect("\\("),
                      gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "",
                           colnames(Xs), perl=TRUE),
                      colnames(Xs))

  terms <- unique(terms_vec)
  terms <- terms[-which(terms=="Intercept")]

  allterminds <- lapply(terms, function(x) which(x == terms_vec))

  ## Loop through terms 
  dfList <- vector("list", length = length(allterminds)) 
  ## gammaSElist <- vector("list", length = length(allterminds)) 

  gammaTermsList <- vector("list", length = length(allterminds))

  for (k in 1:length(allterminds)) {

    terminds <- allterminds[[k]]

    Xs_sub <- matrix(Xs[, terminds], ncol = length(terminds))
    Q_sub <- matrix(Q[terminds, ], nrow = length(terminds))

    

    gammaTermSamples <- Xs_sub %*% Q_sub

    gammaTermsList[[k]] <- gammaTermSamples

    ## gammaTermSamples <- Xs[, terminds] %*% Q[terminds, ]

    ## if (length(terminds )==1) {
    ##   Xs_sub_int <- matrix(Xs[, 1], ncol = 1)
    ##   Q_sub_int <- matrix(Q[1, ], nrow = 1)
    ##   gammaTermSamples <- gammaTermSamples - Xs_sub_int %*% Q_sub_int
    ## }

    gammaTermMean <- rowMeans(gammaTermSamples)
    gammaTermLo <- apply(gammaTermSamples, 1, quantile, probs = alpha / 2)
    gammaTermHi <- apply(gammaTermSamples, 1, quantile, probs = 1 - alpha / 2)



    ## dfList[[k]] <- data.frame(
    ##   j = terms[k],
    ##   xj = df %>% pull(terms[k]),
    ##   fjxj = gammaTermMean,
    ##   lo = gammaTermLo,
    ##   hi = gammaTermHi
    ## )

    dfList[[k]] <- data.frame(
      j = terms[k],
      xj = df %>% pull(terms[k]),
      fjxj = gammaTermMean - mean(gammaTermMean),
      lo = gammaTermLo - mean(gammaTermMean),
      hi = gammaTermHi - mean(gammaTermMean)
    )

    ## gammaSElist[[k]] <- data.frame(
    ##   j = terms[k],
    ##   xj = x[, k],
    ##   fjxj = gammaTerms[[1]][, k],
    ##   lo = gammaTerms[[1]][, k] - 1.96 * gammaTerms[[2]][, k],
    ##   hi = gammaTerms[[1]][, k] + 1.96 * gammaTerms[[2]][, k]
    ## )
    
    ## print(k)

    
  }

  names(gammaTermsList) <- terms

  gammaDf <- dfList %>% plyr::rbind.fill()

  ## Make a plot
  gamPlot <- gammaDf %>%
    ggplot() +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_rug(aes(xj, fjxj), sides = "b", alpha = 0.1) +
    ## geom_line(aes(xj, lo), col = "firebrick3", lty = "dashed") +
    ## geom_line(aes(xj, hi), col = "firebrick3", lty = "dashed") +
    geom_ribbon(aes(xj, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.5) +
    geom_line(aes(xj, fjxj), col = "firebrick3") +
    facet_wrap(~j, scales = "free_x") +
    labs(title = "Projected additive summary of GP fit",
         subtitle = "Using posterior draws of GP",
         x = "",
         y = TeX("$h_j(x_j)$")) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 

  ## Output
  list(gamFit = gamFit,
       gamDf = gammaDf,
       gamPlot = gamPlot,
       gammaTerms = gammaTermsList,
       call = gamcall)


  
}
