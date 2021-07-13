
##' Plots the results of an additive summary
##'
##' Plots the partial/adjusted effect curves from an additive summary
##' @title Plot an additive summary
##' @param additive_summary
##' @param ribbonFill 
##' @return
##' @author Spencer Woody
additive_summary_plot <- function(additive_summary, ribbonFill = "grey80") {

  additive_summary$gamDf %>%
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
