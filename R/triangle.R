##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Triangles
##' @param x
##' @param fsamples
##' @param num_quantiles
##' @return
##' @author Spencer Woody
triangle <- function(x, fsamples, num_quantiles = 100) {
  
  if (length(unique(x)) < num_quantiles) {
    warning(sprintf("num_quantiles argument (%i) less than unique values of x (%i)", num_quantiles, length(unique(x))))
    num_quantiles <- length(unique(x))
    #  myinds <- order(x)
    temp_df2 <- data.frame(x = x) %>% mutate(Row = row_number()) %>% arrange(x) %>% distinct(x, .keep_all=TRUE)
  } else {
    temp_df <- data.frame(x = x) %>% mutate(Row = row_number())
    
    myquants <- quantile(x, probs = seq(0.025, 0.975, length.out = num_quantiles), type = 3)
    
    temp_df2 <- temp_df %>% filter(x %in% myquants) %>% arrange(x) %>% distinct(x, .keep_all = TRUE)
    
    if (nrow(temp_df2) < num_quantiles)  {
      warning("Fewer quantiles reported due to ties")
      num_quantiles <- nrow(temp_df2)
    }
    
    
    myinds <- which(x %in% myquants)
  }
  
  myrowvec <- temp_df2 %>% pull(Row)
  xvec <- temp_df2 %>% pull(x)
  
  ## print(num_quantiles)
  ## print(myrowvec)
  ## print(max(myrowvec))
  ## print(dim(fsamples))
  
  veclen <- num_quantiles * (num_quantiles - 1) / 2
  
  ## print(veclen)
  
  l_vec <- rep(NA, veclen)
  u_vec <- rep(NA, veclen)
  x_l_vec <- rep(NA, veclen)
  x_u_vec <- rep(NA, veclen)
  diff_vec <- rep(NA, veclen)
  prob_vec <- rep(NA, veclen)
  
  idx <- 1
  for (l in 1:(num_quantiles - 1)) {
    for (u in (l + 1):num_quantiles) {
      l_vec[idx] <- l/num_quantiles
      u_vec[idx] <- u/num_quantiles
      x_l_vec[idx] <- xvec[l][1]
      x_u_vec[idx] <- xvec[u][1]
      diff_vec[idx] <- mean(fsamples[myrowvec[u], ] - fsamples[myrowvec[l], ])
      prob_vec[idx] <-  mean(fsamples[myrowvec[u], ] > fsamples[myrowvec[l], ])
      idx <- idx+1
    }
  }
  
  data.frame(
    l = l_vec,
    u =u_vec,
    x_l =x_l_vec,
    x_u = x_u_vec,
    diff = diff_vec,
    prob = prob_vec)
  
}



