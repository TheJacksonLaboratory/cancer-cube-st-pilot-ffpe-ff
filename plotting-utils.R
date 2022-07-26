#' Add a title to a cowplot
#' 
#' A convenience function that adds a panel with a title above an existing figure.
#'
#' @param g A ggplot.
#' @param title A string title.
#' @param size Size of title text.
#' @return A ggplot
add.title.to.plot <- function(g, title, size = 14) {
  g.title <- ggdraw() +  draw_label(title, size = size)
  g <- cowplot::plot_grid(title = g.title, g, nrow=2, rel_heights=c(0.1,1))
  g
}

#' Create a plot of the raw H&E image
#' 
#' A convenience function that wraps Seurat's SpatialFeaturePlot to plot only the overlaid H&E image and not the spots.
#'
#' @param obj A Seurat object.
#' @return a ggplot
plot.hne <- function(obj, keep.invisible.legend = FALSE) {
  # Set the opacity/alpha to 0 so that we only see the H&E image.
  g <- SpatialFeaturePlot(obj, features = "nCount_Spatial", alpha = c(0,0))
  # g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
  if(keep.invisible.legend) {
    g <- g + theme(legend.position = "bottom") 
    g <- g + theme(legend.title = element_text(color = "transparent"), legend.text = element_text(color = "transparent"))
    # g <- g + scale_fill_continuous(fill = guide_legend(override.aes = list(alpha = 1)))
    # g <- g + scale_fill_gradientn(colors = brewer.pal(10,"Spectral"), guide = guide_legend(override.aes = list(alpha = 1)))
    # g <- g + guides(fill = guide_legend(override.aes = list(alpha = 1)))
    # See https://stackoverflow.com/questions/45998396/unset-existing-scale-fill-discrete-in-ggplot2-or-suppress-message-for-new-scale
    # for the following code, which suppresses 
    # "Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale."
    i <- which(sapply(g$scales$scales, function(x) 'fill' %in% x$aesthetics))
    g$scales$scales[[i]] <- NULL
    g <- g + scale_fill_gradientn(colours = c("white"))
  } else {
    g <- g + theme(legend.position = "none")
  }
  g
}

#' Create a panel of plots, each showing deconvolved fractions of a particular population.
#' 
#' @param rctd An RCTD object, from the spacexr package.
#' @param pop.weights A data frame whose rows are spots, whose columns are deconvolved populations, and whose 
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    Such a data.frame can be created by format.rctd.output_.
#' @param populations A vector of populations (subset of columns of pop.weights) to plot.
#' @param obj A SeuratObject holding the spatial features data. Only required if show.hne = TRUE
#' @param show.hne Boolean indicating whether to also plot the H&E.
#' @param use.absolute.scale Boolean indicating whether the legend of each subplot should range from 0 to 1.
#' @return A ggplot
plot.population.weights <- function(rctd, pop.weights, populations, obj = NULL, show.hne = TRUE, use.absolute.scale = FALSE, title.size = 20) {
  if(use.absolute.scale) { p_load(ggpubr) }
  names(populations) <- populations
  plts <-
    llply(populations,
          .fun = function(pop) {
            mx <- max(pop.weights[,pop])
            if(use.absolute.scale) { mx <- 1 }
            g <- plot_puck_continuous(rctd@spatialRNA, colnames(rctd@spatialRNA@counts), as.matrix(pop.weights)[, pop], size = 2, ylimit = c(0, mx))
            g <- g + theme_void() + scale_y_reverse()
            g <- g + labs(color="fraction")
            g <- g + ggtitle(pop) + theme(plot.title = element_text(size=title.size))
            g
          })
  g.legend <- get_legend(plts[[1]])
  if(use.absolute.scale) { 
    plts <- llply(plts, .fun = function(g) g + theme(legend.position = "none"))
  }
  ncol = 2
  if(show.hne) {
    g.hne <- plot.hne(obj)
    g.hne <- g.hne + ggtitle("H&E") + theme(plot.title = element_text(size=title.size))
    xrange <- ggplot_build(g.hne)$layout$panel_scales_x[[1]]$range$range
    yrange <- ggplot_build(g.hne)$layout$panel_scales_y[[1]]$range$range
    # g.hne <- g.hne + labs(title = element_blank())
    p_load(egg) # for set_panel_size
    # plts <- lapply(plts, set_panel_size, width = unit(xrange[2] - xrange[1], "npc"), height = unit((yrange[2]- yrange[1]) + 100, "npc"))
    plts <- c(list(g.hne), plts)
    ncol = 3
    g <- plot_grid(plotlist=lapply(plts, set_panel_size, width = unit(xrange[2] - xrange[1], "null"), height = unit(yrange[2]- yrange[1], "null")), ncol=ncol)
    if(use.absolute.scale) {
      g <- plot_grid(g, g.legend, nrow = 1, rel_widths = c(10, 0.75))
    }
    return(g)
  }
  if(use.absolute.scale) {
    g <- plot_grid(plotlist = plts, ncol=ncol)
    g <- plot_grid(g, g.legend, nrow = 1, rel_widths = c(10, 0.75))
  }
  return(g)
}
