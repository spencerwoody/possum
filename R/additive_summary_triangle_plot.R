##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Triangle
##' @return
##' @author Spencer Woody
##'
##' @export
additive_summary_triangle_plot <- function(possum1, windsor=NA) {

  temp<-possum1$triangleDf
  
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
