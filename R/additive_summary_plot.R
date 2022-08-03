
##' Plot additive summary
##'
##' Plots the partial/adjusted effect curves from additive_summary
##' @title Plot an additive summary
##' @param additive_summary Output of additive_summary
##' @param ribbonFill An optional argument for the color of the credible interval
##' @param windsor "Windsorized" summary, which removes the top and bottom windsor/2 quantiles of each covariate. Only available if "quants" present in the posterior summary
##' @return
##' @export
##' @author Spencer Woody
additive_summary_plot <- function(additive_summary, ribbonFill = "grey80",
                                  windsor=NA) {

  temp<-additive_summary$gamDf
  
  if (!is.na(windsor)) {
    if (!("quant" %in% colnames(temp))) {
      stop("Quantiles not supplied") 
    }
    temp <- temp %>% filter(quant > windsor/2 & quant < 1-windsor/2)
    glimpse(temp)
  }
  
   temp %>%
    distinct() %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(x_j, ymin = fx_j_lo, ymax = fx_j_hi),
                fill = ribbonFill, alpha = 0.5) +
    geom_line(aes(x_j, fx_j_mean), col = "firebrick3") +
    geom_rug(aes(x_j, fx_j_mean), sides = "b", alpha = 0.25) +
    facet_wrap(~term, scale = "free_x") +
                                        # theme_minimal(base_size = 14) +
                                        # theme_half_open() +
    ## theme_cowplot() +
                                        # labs(x = TeX("$x_j$"), y = TeX("$f_j(x_j)$")) +
    labs(x = ("term"), y = ("Partial effect"))
}
