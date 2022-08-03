##' Plot summary gradient
##'
##' Create plot of posterior probability for positive/negative gradient in the summary curves
##'
##' @param additive_summary Output of additive_summary
##' @param windsor windsor "Windsorized" summary, which removes the top and bottom windsor/2 quantiles of each covariate
##'
##' @title Plot summary gradient
##' @return
##' @author Spencer Woody
##'
##' @export
additive_summary_triangle_plot <- function(additive_summary, windsor=NA) {

  temp<-additive_summary$triangleDf
  
  if (!is.na(windsor)) {
    temp <- temp %>% 
      filter(u > windsor/2 & u < 1-windsor/2) %>% 
      filter(l > windsor/2 & l < 1-windsor/2)
  }
  
  temp %>%
  ggplot() +
  geom_raster(aes(l, u, fill = prob)) +
  facet_wrap(~term) +
  scale_fill_gradientn(limits = c(0,1),
                         colors = brewer.pal(7, "RdBu"),
                         values = c(1, 0.85, 0.75, 0.5, 0.25, 0.15, 0)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=0.5)) +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5)) +
  coord_equal() +
  labs(x = "Lower quantile", y = "Upper quantile")

}
