##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Triangle
##' @return
##' @author Spencer Woody
##'
##' @export
additive_summary_triangle_plot <- function(possum1) {

  possum1$triangleDf %>%
  ggplot() +
  geom_raster(aes(l, u, fill = prob)) +
  facet_wrap(~term) +
  scale_fill_gradient2("", midpoint=0.5) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=0.5)) +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5)) +
  coord_equal() +
  labs(x = "Lower quantile", y = "Upper quantile")

}
