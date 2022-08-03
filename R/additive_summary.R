
make_grid <- function(df, quants=seq(0, 1, by=0.005)) {
  require(dplyr)
  
  pen_num <- df %>% select(where(is.numeric))
  pen_fac <- df %>% select(where(is.factor))
  
  pen_num_levels <- apply(pen_num, 2, quantile, probs=quants) %>% as.data.frame()
  
  if (ncol(pen_fac)==0) {
    pen_levels <- pen_num_levels %>% select(colnames(df))
  } else {
    pen_fac_levels <- dplyr::distinct(pen_fac) %>% as.data.frame()
    n1 <- nrow(pen_num_levels)
    n2 <- nrow(pen_fac_levels)
    pen_fac_levels <- pen_fac_levels %>% 
      slice(rep(1:n(), each = ceiling(n1/n2))) %>% 
      slice(1:n1)
    pen_levels <- cbind(pen_num_levels, pen_fac_levels) %>% select(colnames(df))
  }
  
  pen_levels
  
}

##' Compute posterior additive summary
##'
##' This function computes the point estimate and credible intervals for the summary of the function f. At a minimum, the user must specify the form of the summary, a matrix of posterior draws of f, and the dataframe which contains the inputs of f
##' @title additive_summary
##'
##' @param summaryCall A gam fomula for the additive summary to be
##'   computed.  Should be in the form of fhat ~ s(x1) + s(x2) + ....
##'   See ?mgcv::formula.gam
##' @param fhatSamples N \times NMC matrix of posterior draws of the function f, where N is the number of observations and NMC is the number of Monte Carlo posterior samples 
##' @param fhat A point estimate (posterior mean) for the function f
##' @param df The dataframe from which the summary will be computed. This should include all the inputs of f
##' @param grid_size 
##' @param alpha The function will return the alpha/2 and 1-alpha/2 credible intervals for the summary. 
##' @param fast If TRUE, the function will compute the summary on a grid, specifically for quantiles of the covariates as specified in quants, of the data rather than for the whole dataset
##' @param quants The quantiles of the covariates in df on which to compute the summary when fast=TRUE
##' @param verbose If TRUE, the function will print out the progress of summary computation
##' @param return_samples If TRUE, the function will of the design matrix for the summary 
##' @param meta A tag for the dataframe of the summary 
##'
##' @return
##' @author Spencer Woody
##' @export
additive_summary <- function(summaryCall,
                             fhatSamples,
                             fhat = rowMeans(fhatSamples),
                             df = NULL,
                             ## terms = NULL,
                             alpha = 0.05,
                             fast=TRUE, quants=seq(0, 1, by=0.005),
                             grid_size = 100,
                             verbose=FALSE,
                             return_samples = TRUE, 
                             meta = NA
)## function(fhatmat, df, gamFit, y, meta=NA,
## alpha = 0.10, ribbonFill = "grey50")
{
  
  if (fast) {
    df_est <- make_grid(df, quants=quants)
  } else {
    df_est <- df
  }
  
  
  ## TODO: improve formula interface by using update: 
  #update(. ~ u+v, res  ~ . ) #> res ~ u + v
  
  
  if (sapply(df, is.factor) %>% any()) {
    warning("Looks like there are factors in your data.frame.\n 
            If so, please be sure to name them as factors explicitly in summary call")
  }
  
  ## Rename... (y not needed?)
  
  if (verbose) cat("Calculating point estimate for the summary...\n")
  
  gamFit <- gam(summaryCall, data = df)
  
  
  ## Regular
  
  if (verbose) cat("Extracting terms of summary...\n")
  
  gamFitTerms <- predict(gamFit, type="terms", newdata=df_est)
  
  Xgam <- model.matrix(gamFit, newdata=df_est)
  
  Xgam_full <- model.matrix(gamFit)
  V <- vcov(gamFit, dispersion = 1)
  Q <- crossprod(V, crossprod(Xgam_full, fhatSamples))
  
  ## S <- Xgam %*% crossprod(V, t(Xgam))
  
  fittedValues <- Xgam_full %*% Q
  
  gamTerms <- attr(gamFit$terms, "term.labels")
  
  gamFitMat <- Xgam %*% Q
  
  ## Loop through all the terms and calculate the credible bands
  gamDfList <- list()
  gamTerm <- list()
  gamFactorDfList <- list()
  
  ## Loop through factor
  
  if (verbose) cat("Looping through factors...\n")
  
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
  
  if (verbose) cat("Looping through all other terms...\n")
  
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
      x_j = df_est[, gamTerms[j]],
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
      # fx_res=fhat-(rowSums(gamFitTerms[, !str_detect(colnames(gamFitTerms), 
      #                                                gamTerms[j]),drop=FALSE]) + 
      #                coef(gamFit)["(Intercept)"]),
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
    myx <- df_est %>% dplyr::pull(triangle_terms[jay])
    triangleDfList[[jay]] <- triangle(myx, gamTerm[[triangle_terms[jay]]]) %>% 
      mutate(term = triangle_terms[jay])
  }
  
  triangleDf <- triangleDfList %>% plyr::rbind.fill()
  
  # cat("Correlation of summary with original fit:", round(sqrt(summary(gamFit)$r.sq, 3), "\n"))
  
  cat(dim(fittedValues))
  cat(dim(fhatSamples))
  
  ## Output
  out = list(
    gamFit = gamFit,
    gamFitMat = gamFitMat, #
    gamDf = gamDf,
    ## gamPlot = gamPlot,
    gamFactorDf = gamFactorDf,
    fittedValues = fittedValues,
    summaryRsq = rsqGamma(fittedValues, fhatSamples),
    Xgam = Xgam,
    Q = Q,
    gamFitTerms = gamFitTerms,
    triangleDf = triangleDf,
    df=df,
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

