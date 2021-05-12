
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot
##' @param possum1
##' @return
##' @author Spencer Woody
possum_plot1 <- function(possum1) {
  ribbonFill <- "grey80"

  possum1$gamDf %>%
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
