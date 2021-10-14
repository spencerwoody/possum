##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title additive_summary
##' @param summaryCall A gam fomula for the additive summary to be
##'   computed.  Should be in the form of fhat ~ s(x1) + s(x2) + ....
##'   See ?mgcv::formula.gam
##' @param fhatSamples
##' @param fhat
##' @param df
##' @param ribbonFill
##' @param meta
##' @param verbose
##' @return
##' @author Spencer Woody
additive_summary_local <- function(summaryCall,
                                   fhatSamples,
                                   fhat = rowMeans(fhatSamples),
                                   df = NULL,
                                   ## terms = NULL,
                                   grid_size = 100,
                                   alpha = 0.05,
                                   ribbonFill = "grey50",
                                   return_samples = TRUE,
                                   meta = NA,
                                   verbose=FALSE
)## function(fhatmat, df, gamFit, y, meta=NA,
## alpha = 0.10, ribbonFill = "grey50")
{
  
  ## TODO: improve formula interface by using update: update(. ~ u+v, res  ~ . ) #> res ~ u + v
  
  
  ## if (sapply(df, is.factor) %>% any()) warning("Looks like there are factors in your data.frame.\n If these, please be sure to explicitly ")
  
  ## Rename... (y not needed?)
  
  if (verbose) cat("Calculating point estimate for the summary...\n")
  
  gamFit <- gam(summaryCall, data = df)
  
  
  ## Regular
  
  if (verbose) cat("Extracting terms of summary...\n")
  
  gamFitTerms <- predict(gamFit, type="terms")
  
  Xgam <- model.matrix(gamFit)
  V <- vcov(gamFit, dispersion = 1)
  Q <- crossprod(V, crossprod(Xgam, fhatSamples))
  
  ## S <- Xgam %*% crossprod(V, t(Xgam))
  
  fittedValues <- Xgam %*% Q
  
  gamTerms <- attr(gamFit$terms, "term.labels")
  
  gamFitMat <- Xgam %*% Q
  
  ## Loop through all the terms and calculate the credible bands
  gamDfList <- list()
  gamTerm <- list()
  gamFactorDfList <- list()
  
  ## Loop through factor
  
  cat("Looping through factors...\n")
  
  if (any(str_detect(gamTerms, "factor"))) {
    gamFactorTerms <- gamTerms[str_detect(gamTerms, "factor")]
    
    for (k in 1:length(gamFactorTerms)) {
      
      termjidx <- which(str_detect(colnames(Xgam), fixed(gamFactorTerms[k])))
      
      termnameK <- str_extract(gamFactorTerms[k], "\\([^()]+\\)") %>%
        str_sub(start = 2, end = -2)
      
      colnamesK <- colnames(Xgam)[termjidx] %>%
        str_remove(fixed(gamFactorTerms[k]))
      
      termKdfList <- vector("list", length(termjidx))
      
      for (ell in 1:length(termjidx)) {
        termKdfList[[ell]] <- data.frame(term = termnameK,
                                         level = colnamesK[ell],
                                         post = Q[termjidx[ell], ],
                                         stringsAsFactors = FALSE)
      }
      
      termKdf <- plyr::rbind.fill(termKdfList)
      
      gamFactorDfList[[k]] <- termKdf
      
      
    }
    
  }
  gamFactorDf <- gamFactorDfList %>% plyr::rbind.fill()
  
  ## Loop through all other terms
  
  cat("Looping through all other terms...\n")
  
  for (j in 1:length(gamTerms)) {
    
    if (verbose) cat(sprintf("Term %i out of %i...\n", j, length(gamTerms)))
    
    if (str_detect(gamTerms[j], "factor")) {
      next
    }
    
    termjidx <- which(str_detect(colnames(Xgam), gamTerms[j]))
    
    # TODO: finish grid spec
    # gamdatix = which(str_detect(colnames(gamFit$model), gamTerms[j]))
    # xx = possum_fit$gamFit$model[,gamdatix]
    # ngr = grid_size#min(grid_size, length(unique(xx)))
    # Xgr = do.call("rbind", replicate(ngr, gamFit$model[1,], simplify = FALSE))
    # Xgr[,gamdatix] = .cp_quantile(xx, ngr)
    
    
    ## Check if term is linear
    if (length(termjidx) == 1) {
      gamTermjPost <- matrix(Xgam[, termjidx], ncol = 1) %*%
        matrix(Q[termjidx, ], nrow = 1)
      
      ## Set so that the Mean is 0...
      gamTermjPost <- gamTermjPost - mean(gamTermjPost)
      
    } else {
      gamTermjPost <- Xgam[, termjidx] %*% Q[termjidx, ]
    }
    
    gamTerm[[j]] <- gamTermjPost
    
    # gamDfLongList[[length(gamDfList)]]
    
    gamDfList[[length(gamDfList) + 1]] <- data.frame(
      term = gamTerms[j],
      x_j = df[, gamTerms[j]],
      ## fx_j_mean = gamFitTerms[, which(str_detect(colnames(gamFitTerms), gamTerms[j]))],
      fx_j_mean = rowMeans(gamTermjPost),
      fx_j_lo = apply(gamTermjPost, 1, quantile, probs = alpha / 2),
      fx_j_hi = apply(gamTermjPost, 1, quantile, probs = 1 - alpha / 2),
      fx_j_lo95 = apply(gamTermjPost, 1, quantile, probs = 0.05 / 2),
      fx_j_hi95 = apply(gamTermjPost, 1, quantile, probs = 1 - 0.05 / 2),
      fx_j_lo80 = apply(gamTermjPost, 1, quantile, probs = 0.20 / 2),
      fx_j_hi80 = apply(gamTermjPost, 1, quantile, probs = 1 - 0.20 / 2),
      fx_j_lo50 = apply(gamTermjPost, 1, quantile, probs = 0.50 / 2),
      fx_j_hi50 = apply(gamTermjPost, 1, quantile, probs = 1 - 0.50 / 2),
      fx_res=fhat-(rowSums(gamFitTerms[, !str_detect(colnames(gamFitTerms), gamTerms[j]),drop=FALSE]) + coef(gamFit)["(Intercept)"]),
      meta = meta
    )
    
  }
  
  names(gamTerm) <- gamTerms
  
  if (verbose) cat("Combining summary components...\n")
  
  gamDf <- plyr::rbind.fill(gamDfList)
  
  
  
  ## Triangle
  triangle_terms <- gamTerms
  triangle_terms <- triangle_terms[!str_detect(triangle_terms, "factor")]
  
  num_triangles <- length(triangle_terms)
  
  triangleDfList <- vector("list", length=num_triangles)
  
  if (verbose) cat("Making triangle heatmaps...\n")
  
  for (jay in 1:num_triangles) {
    if (verbose) cat (sprintf("Triangle %i out of %i...\n", jay, num_triangles))
    myx <- df %>% dplyr::pull(triangle_terms[jay])
    triangleDfList[[jay]] <- triangle(myx, gamTerm[[triangle_terms[jay]]]) %>% mutate(term = triangle_terms[jay])
  }
  
  triangleDf <- triangleDfList %>% plyr::rbind.fill()
  
  # cat("Correlation of summary with original fit:", round(sqrt(summary(gamFit)$r.sq, 3), "\n"))
  
  ## Output
  out = list(
    gamFit = gamFit,
    gamFitMat = gamFitMat, #
    gamDf = gamDf,
    ## gamPlot = gamPlot,
    gamFactorDf = gamFactorDf,
    fittedValues = fittedValues,
    summaryRsq = rsqGamma(gamFitMat, fhatSamples),
    Xgam = Xgam,
    Q = Q,
    gamFitTerms = gamFitTerms,
    triangleDf = triangleDf,
    call = summaryCall
  )
  
  if(return_samples) out$gamTermPostSamples = gamTerm
  
  return(out)
  
}


.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))
  
  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/(nobs-1)))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }
  
  return(ret)
}

