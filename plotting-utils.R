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

#' Plot the distribution of fractions for each (deconvolved population) across samples.
#' 
#' If there is only one sample, populations will be listed on the x axis. 
#' Otherwise, plots will be faceted on populations and samples will be listed on x axis.
#' 
#' @param pop.weights A data frame whose rows are spots, whose columns are deconvolved populations, and whose 
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    data.frames for each sample can be created by format.rctd.output_.
#' @param sample.col The column within pop.weights that gives the name of the sample.
#'                   If null, pop.weights is considered to describe only one sample.
#' @return a ggplot
plot.population.fractions.across.samples <- function(pop.weights, sample.col = NULL) {
  num.samples <- 1
  if(!is.null(sample.col) && (length(unique(pop.weights[, sample.col])) > 1)) {
    df <- reshape2::melt(pop.weights, id.vars = sample.col)
    colnames(df) <- c(sample.col, "variable", "value")
    df <- subset(df, !(variable %in% c("x","y")))
    g <- ggplot(data = df, aes_string(x = sample.col, y = "value"))
    #g <- g+scale_fill_manual(values=colSide)
    g <- g + facet_wrap(as.formula(paste("~", "variable")), scales = "free_y")
    g <- g + xlab("Sample")
  } else {
    df <- reshape2::melt(pop.weights)
    colnames(df) <- c("variable", "value")
    df <- subset(df, !(variable %in% c("x","y")))
    g <- ggplot(data = df, aes(x = variable, y = value))
    g <- g + xlab("Population")
  }
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  # beeswarm takes a long time
  # g <- g + geom_beeswarm()
  # g <- g + geom_boxplot()
  g <- g + geom_violin(trim=FALSE)
  g <- g +scale_fill_brewer(palette="Dark2")
  #g +scale_fill_manual(values=c("#d8b365", "#d8b365", "#5ab4ac", "#5ab4ac"))
  g <- g + ylab("Population Fraction")
  g
}

#' Plot distribution of a response variable (e.g., UMI count) across spots relative to a dependent variable (e.g., tissue vs background spot status)
#' 
#' Create violin/boxplots showing distribution across spots of a user-specified response variable relative to a dependent variable and facted by sample.
#'
#' @param df A data.frame holding the response, depenndent, and faceting variables.
#' @param response.var The variable whose distribution over spots will be plotted.
#' @param facet.var The faceting variable.
#' @param dependent.var The variable that the response variable will be plotted as a function of.
#' @return A ggplot
plot.distributions.vs.cell.type <- function(df, response.var = "nCount_Spatial", facet.var = "sample.name", dependent.var = "spot_type") {
  lvls <- unique(df[, facet.var, drop=TRUE])
  df[,facet.var] <- factor(df[,facet.var], levels=lvls)
  stat.test <- df %>% group_by_at(facet.var) %>% wilcox_test(as.formula(paste0(response.var, " ~ ", dependent.var)), p.adjust.method = "none")
  stat.test$p.adj <- p.adjust(stat.test$p, method = "bonferroni")
  stat.test <- stat.test %>% add_xy_position(x = dependent.var)
  stat.test$y.position <- log2(stat.test$y.position)
  stat.test$label <- stars.pval(stat.test$p.adj)
  g <- ggviolin(df, dependent.var, response.var, facet.by = facet.var, add = "boxplot")
  g <- g + yscale("log2") + stat_pvalue_manual(stat.test, label="label") + xlab("Spot Type")
  g
}

#' Create a spatial feature plot
#' 
#' A convenience function that wraps Seurat's SpatialFeaturePlot to add a centered title and a legend at the bottom.
#'
#' @param obj A Seurat object.
#' @param features A vector of strings listing one or more features to plot.
#' @param legend.name A string name for the legend
#' @param rescale.legend Boolean indicating whether to divide the values by 1000 in legend scale (obviously hacky, but intended for large numbers like total counts)
#' @return a ggplot
plot.spatial <- function(obj, features = c("nCount_Spatial"), legend.name = "Read Count", rescale.legend = TRUE) {
  g <- SpatialFeaturePlot(obj, features = features)
  # g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15), legend.position = "bottom", legend.key.width = unit(1.5, 'cm'))
  g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15), legend.position = "bottom")
  # See https://stackoverflow.com/questions/45998396/unset-existing-scale-fill-discrete-in-ggplot2-or-suppress-message-for-new-scale
  # for the following code, which suppresses 
  # "Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale."
  # g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0fk', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  if(rescale.legend) {
    i <- which(sapply(g$scales$scales, function(x) 'fill' %in% x$aesthetics))
    g$scales$scales[[i]] <- NULL
    g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0f', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  }
  g
}

#' Create a spatial feature plot
#' 
#' A convenience function that wraps Seurat's SpatialFeaturePlot to add a centered title and a legend at the bottom.
#'
#' @param obj A Seurat object.
#' @param features A vector of strings listing one or more features to plot.
#' @param include.hne Boolean indicating whether to include an H&E plot
#' @param include.umi.cnts Boolean indicating whether to include a plot of UMI counts
#' @param include.feature.cnts Boolean indicating whether to include a plot of feature / gene counts
#' @return a ggplot
plot.features <- function(obj, features, include.hne = FALSE, include.umi.cnts = FALSE, include.feature.cnts = FALSE) {
  names(features) <- features
  plts <- lapply(features, function(feature) plot.spatial(obj, features = c(feature), legend.name = feature, rescale.legend = FALSE))
  if(include.feature.cnts) {
    p <- plot.spatial(obj, "nFeature_Spatial", "# Features (K)", rescale.legend = TRUE)
    plts <- c(list(p), plts)
  }
  if(include.umi.cnts) {
    p <- plot.spatial(obj, "nCount_Spatial", "# UMIs (K)", rescale.legend = TRUE)
    plts <- c(list(p), plts)
  }
  if(include.hne) {
    p <- plot.hne(obj, keep.invisible.legend = TRUE)
    plts <- c(list(p), plts)
  }
  plot_grid(plotlist = plts)
}

#' Create a VennDiagram plot from a list on genesets
#' @param temp_list list with the genesets (1-4) we want to plot in the VennDiagram
#' @param title Title of the plot
#' @return a ggplot
plot.VennDiagram.list <- function(temp_list, title){
  g <- ggVennDiagram(temp_list)
  g <- g + ggtitle(title) 
  g <- g + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  g <- g + scale_x_continuous(expand = expansion(mult = .2))
  g
}

#' Perform GO and KEGG analysis using clusterProfiler library
#' 
#' @param geneset A list with genes
#' @param analysis_set Prefix character describing the geneset
GO.KEGG.enirhcment.analysis<-function(geneset,analysis_set){
  yy <- enrichGO(geneset, OrgDb=org.Hs.eg.db, keyType= 'SYMBOL', ont = "ALL", pvalueCutoff=0.01, qvalueCutoff = 0.05)
  png(paste0(plots_dir, "/", analysis_file_prefix, "GO_gene_analysis_dotplot_", analysis_set ,".png"), width = 800, height = 1600)
  dotplot(yy, split="ONTOLOGY",color = "qvalue",showCategory = 15) + facet_grid(ONTOLOGY~., scale="free")
  d <- dev.off()
  
  querry_str="detec"
  yy@result[["Description"]][grepl(querry_str , yy@result[["Description"]], fixed = TRUE)]
  
  ontol_set=c("BP","MF","CC")
  for (ontol in ontol_set){
    yy <- enrichGO(geneset, OrgDb=org.Hs.eg.db, keyType= 'SYMBOL', ont = ontol, pvalueCutoff=0.01, qvalueCutoff = 0.05)
    png(paste0(plots_dir, "/", analysis_file_prefix, "GO_gene_analysis_goplot_",ontol,"_" ,analysis_set ,".png"), width = 800, height = 1600)
    goplot(yy)
    d <- dev.off()
  }
  
  gene.df <- bitr(geneset, fromType = "SYMBOL", toType = c("ENTREZID" ), OrgDb = org.Hs.eg.db)
  mkk <- enrichKEGG(gene = gene.df$ENTREZID , organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1)
  png(paste0(plots_dir, "/", analysis_file_prefix, "KEGG_dotplot_",ontol,"_" ,analysis_set ,".png"), width = 800, height = 1600)
  dotplot(mkk, title = analysis_set)
  d <- dev.off()
}

#' Create a panel of plots, each of which shows the value of a feature (gene or metadata) on the y axis
#' with the spots linearized on the x axis.
#' 
#' @param obj A Seurat object.
#' @param features A vector of strings listing one or more features to plot.
#'                 Each feature should be a row in the "Spatial" assay or a column in the object's metadata.
#' @param order.by A vector of strings listing a subset of features used to order the spots on the x axis.
#' @return a ggplot
create.feature.strip.plot <- function(obj, features, order.by = NULL) {
  mat <- Seurat::GetAssayData(obj, assay="Spatial")
  expr <- cpm(as.matrix(mat), log = FALSE)
  meta <- obj[[]]
  common.spots <- intersect(rownames(meta),colnames(mat))
  mat <- mat[, common.spots]
  meta <- meta[common.spots, ]
  merged <- rbind(as.data.frame(mat[rownames(mat) %in% features,]), t(meta[, colnames(meta) %in% features]))
  m <- reshape2::melt(as.matrix(merged))
  colnames(m) <- c("feature", "sample", "value")
  if(!is.null(order.by)) {
    tmp <- rbind(mat[rownames(mat) %in% order.by,], t(meta[, colnames(meta) %in% order.by]))
    order.by <- order.by[order.by %in% rownames(tmp)]
    tmp <- tmp[order.by,]
    ii <- do.call('order', as.data.frame(t(tmp)))
    lvls <- colnames(tmp)[ii] 
    m$sample <- factor(m$sample, levels = lvls)
    m$feature <- factor(m$feature, levels = rev(c(order.by, features[!(features %in% order.by)])))
  }
  g <- ggplot(data = m) + geom_col(aes(x = sample, y = value)) + facet_wrap(~feature, scales="free", ncol = 1)
  g <- g + ylab("") + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  return(list("merged" = merged, "g" = g))
}
