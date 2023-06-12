#' Add a title to a cowplot
#' 
#' A convenience function that adds a panel with a title above an existing figure.
#'
#' @param g A ggplot.
#' @param title A string title.
#' @param size Size of title text.
#' @return A ggplot
add.title.to.plot <- function(g, title, size = 14, rel_heights=c(0,1)) {
  g.title <- ggdraw() +  draw_label(title, size = size)
  g <- cowplot::plot_grid(title = g.title, g, nrow=2, rel_heights=rel_heights)
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
    g <- g + theme(legend.position = "right") 
    ## g <- g + theme(legend.title = element_text(color = "transparent"), legend.text = element_text(color = "transparent"))
    ### See https://stackoverflow.com/questions/45998396/unset-existing-scale-fill-discrete-in-ggplot2-or-suppress-message-for-new-scale
    ### for the following code, which suppresses 
    ### "Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale."
    ##i <- which(sapply(g$scales$scales, function(x) 'fill' %in% x$aesthetics))
    ##g$scales$scales[[i]] <- NULL
    ##g <- g + scale_fill_gradientn(colours = c("white"))
    g <- plot_grid(g + theme(legend.position = "none"), NULL, nrow = 1, rel_widths = c(7,3))
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
#' @param id.col The column(s) within pop.weights _other_ than the populations.
#' @param sample.col The column within pop.weights that gives the name of the sample.
#'                   If null, pop.weights is considered to describe only one sample.
#' @return a ggplot
plot.population.fractions.across.samples <- function(pop.weights, id.cols, sample.col = NULL, ...) {
  num.samples <- 1
  if(!is.null(sample.col) && (length(unique(pop.weights[, sample.col])) > 1)) {
    df <- reshape2::melt(pop.weights, id.vars = id.cols)
    colnames(df) <- c(id.cols, "variable", "value")
    df <- subset(df, !(variable %in% c("x","y")))
    g <- ggplot(data = df, aes_string(x = sample.col, y = "value", ...))
    #g <- g+scale_fill_manual(values=colSide)
    g <- g + facet_wrap(as.formula(paste("~", "variable")), scales = "free_y")
    g <- g + xlab("Sample")
  } else {
    df <- reshape2::melt(pop.weights)
    colnames(df) <- c("variable", "value")
    df <- subset(df, !(variable %in% c("x","y")))
    g <- ggplot(data = df, aes(x = variable, y = value, ...))
    g <- g + xlab("Population")
  }
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  # beeswarm takes a long time
  # g <- g + geom_beeswarm()
  g <- g + geom_boxplot()
  # g <- g + geom_violin(trim=FALSE)
  #g <- g +scale_fill_brewer(palette="Dark2")
  #g +scale_fill_manual(values=c("#d8b365", "#d8b365", "#5ab4ac", "#5ab4ac"))
  g <- g + ylab("Population Fraction")
  g
}

#' Plot distribution of a response variable (e.g., UMI count) across spots relative to a dependent variable (e.g., tissue vs background spot status)
#' 
#' Create violin/boxplots showing distribution across spots of a user-specified response variable relative to a dependent variable and facted by sample.
#'
#' @param df A data.frame holding the response, dependent, and faceting variables.
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
plot.spatial <- function(obj, features = c("nCount_Spatial"), legend.name = "Read Count", rescale.legend = TRUE, legend.limits = NULL, slot = "counts") {
  g <- SpatialFeaturePlot(obj, features = features, slot = slot)
  # g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15), legend.position = "bottom", legend.key.width = unit(1.5, 'cm'))
  g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15), legend.position = "right")
  # See https://stackoverflow.com/questions/45998396/unset-existing-scale-fill-discrete-in-ggplot2-or-suppress-message-for-new-scale
  # for the following code, which suppresses 
  # "Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale."
  # g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0fk', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  i <- which(sapply(g$scales$scales, function(x) 'fill' %in% x$aesthetics))
  g$scales$scales[[i]] <- NULL
  if(rescale.legend) {
    g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0f', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  } else if(!is.null(legend.limits)) {
    g <- g + scale_fill_gradientn(name = legend.name, colours = Seurat:::SpatialColors(n = 100), limits = legend.limits)
  } else {
    g <- g + scale_fill_gradientn(name = legend.name, colours = Seurat:::SpatialColors(n = 100))
  }
  # g <- g + theme(legend.text = element_text(angle = 45, vjust = 1, hjust=1))
  g <- g + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  g.leg <- get_legend(g) # + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  plot_grid(g + theme(legend.position = "none"), g.leg, nrow = 1, rel_widths = c(7,3))
  # g
}

#' Create a spatial feature plot
#' 
#' A convenience function that wraps Seurat's SpatialFeaturePlot to add a centered title and a legend at the bottom.
#'
#' @param obj A Seurat object.
#' @param features A vector of strings listing one or more features to plot.
#' @param feature.names A vector of strings to use for the corresponding feature in the legend.
#' @param include.hne Boolean indicating whether to include an H&E plot
#' @param include.umi.cnts Boolean indicating whether to include a plot of UMI counts
#' @param include.feature.cnts Boolean indicating whether to include a plot of feature / gene counts
#' @return a ggplot
plot.features_ <- function(obj, features, feature.names = NULL, slot = "counts", include.hne = FALSE, include.umi.cnts = FALSE, include.feature.cnts = FALSE, legend.limits = NULL, ...) {
  plts <- NULL
  if(!is.null(features) && (length(features) > 0)) {
    if(is.null(feature.names)) { feature.names <- features }
    indices <- 1:length(features)
    names(indices) <- features
    lims <- legend.limits
    if(is.null(legend.limits)) { lims <- rep(NULL, length(features))}
    plts <- lapply(indices, function(i) plot.spatial(obj, features = c(features[i]), slot = slot, legend.name = feature.names[i], rescale.legend = FALSE, legend.limits = lims[[i]]))
  }
  if(include.feature.cnts) {
    p <- plot.spatial(obj, "nFeature_Spatial", "# Genes (K)", rescale.legend = TRUE)
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
  #plot_grid(plotlist = plts, ...)
  plts
}

plot.features <- function(obj, features, feature.names = NULL, slot = "counts", include.hne = FALSE, include.umi.cnts = FALSE, include.feature.cnts = FALSE, legend.limits = NULL, ...) {
  plts <- plot.features_(obj, features, feature.names, slot, include.hne, include.umi.cnts, include.feature.cnts , legend.limits, ...)
  plot_grid(plotlist = plts, ...)
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

#' Create a panel of plots for a given feature, one per sample.
#'  
#' @param objs A named list of Seurat 10X spatial objects, each representing a sample.
#' @param titles A named list of titles for each plot, indexed by the sample name (i.e., name of one of the objs).
#' @param feature A string giving the name of the feature to plot.
#' @param legend.name A string giving the name of the legend in each plot.
#' @return A ggplot
plot.feature.across.samples <- function(objs, titles, feature = "nCount_Spatial", legend.name = "Count") {
  samples <- names(objs)
  names(samples) <- samples
  plts <-   
    llply(samples,
          .fun = function(sample) {
            g <- plot.spatial(objs[[sample]], features = c(feature), legend.name = legend.name)
            g <- g + theme(legend.text=element_text(size=9))
            title <- titles[[sample]]
            g <- g + ggtitle(title)
            g
          })
  g.all <- plot_grid(plotlist = plts)
  g.all
}

plot.gene.expression.distribution.relative.to.quantiles <- function(obj, genes, add.housekeeping = TRUE, probs = seq(0.5, 1, by=0.1), ...) {
  raw.cnt.mat <- GetAssayData(obj, assay = "Spatial", slot = "counts")
  expr.mat <- cpm(raw.cnt.mat, log = TRUE)
  nz.genes <- rowSums(raw.cnt.mat) > 0
  raw.cnt.mat <- raw.cnt.mat[nz.genes, ]
  expr.mat <- expr.mat[nz.genes, ]
  
  # Add housekeeping genes
  if(add.housekeeping) {
    # genes <- unique(c(genes, "ACTB", "GAPDH"))
    genes <- unique(c(genes, "ACTB"))
  }
  
  # qs <- calculate.expression.quantiles(raw.cnt.mat, summary.func = median)
  qs <- calculate.expression.quantiles(expr.mat, summary.func = median)

  # Calculate the quantiles and the genes closest to those quantiles
  gene.summaries <- apply(expr.mat, 1, mean)
  # Let's only plot those above the 50% percentile
  qs <- quantile(gene.summaries, probs = probs)
  qs.names <- unlist(llply(qs, .fun = function(val) names(which.min((gene.summaries-val)^2))[1]))
  qs.names <- as.data.frame(qs.names)
  colnames(qs.names)[1] <- "gene"
  qs.names$quantile <- rownames(qs.names)
  qs.names$label <- paste0(qs.names$gene, " (", qs.names$quantile, ")")
  
  mat <- raw.cnt.mat
  genes <- genes[genes %in% rownames(mat)]
  genes <- names(sort(gene.summaries[genes]))
  #qs.names <- subset(qs.names, !(gene %in% genes))
  genes <- genes[!(genes %in% qs.names$gene)]
  labels <- genes
  genes <- c(genes, qs.names$gene)
  labels <- c(labels, qs.names$label)
  names(labels) <- genes
  
  if(FALSE) {
  summarized.expr = apply(mat, 1, median)
  all.summarized.expr.df <- data.frame(gene = rownames(mat), expr = as.numeric(summarized.expr))
  all.summarized.expr.df <- all.summarized.expr.df[order(all.summarized.expr.df$expr, decreasing=TRUE),]
  g2 <- ggplot(data = all.summarized.expr.df, aes(x = expr))
  g2 <- g2 + geom_density()
  
  gene.df <- reshape2::melt(as.matrix(expr.mat[genes, ]))
  }
  
  gene.df <- reshape2::melt(as.matrix(raw.cnt.mat[genes, ]))
  colnames(gene.df) <- c("gene", "sample", "expr")
  title <- paste0(obj[[]]$orig.ident[1], " (", ncol(mat), " spots)")
  # g1 <- ggplot(data = gene.df, aes(x = expr, y = gene)) + geom_boxplot()
  g1 <- ggplot(data = gene.df, aes(x = expr)) + geom_bar() + facet_wrap(~ gene, scale = "free", labeller = as_labeller(labels), ...)
  g1 <- g1 + ylab("Frequency") + xlab("Gene Read Count (Quantiles in Parentheses)")
  g1 <- g1 + ggtitle(title)
  g1
  
} 

#' Plot the fraction of reads by biotype within each sample.
#'  
#' @param objs A named list of Seurat 10X spatial objects, each representing a sample.
#' @return A ggplot
plot.biotypes.across.samples <- function(objs, species = "human", sample.labels = NULL) {
  if (species == 'human'){
    gene_db = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  } else {
    gene_db = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  frac.df <-
    ldply(objs, .parallel = FALSE,
          .fun = function(obj) {
            mat <- GetAssayData(obj, assay="Spatial", slot="counts")
            df <- as.data.frame(mat)
            df[, "gene"] <- rownames(df)
            gb <- get.biotypes_(df$gene, gene_db)
            gb <- condense.biotypes(gb)
            df <- merge(df, gb, all.x = TRUE)
            cnts <- ddply(df, .variables = c("gene_biotype"),
                          .fun = function(df.biotype) {
                            colSums(df.biotype[, !(colnames(df.biotype) %in% c("gene", "gene_biotype"))])
                          })
            rownames(cnts) <- cnts$gene_biotype
            numeric.cols <- colnames(cnts)[!(colnames(cnts) %in% c("gene_biotype"))]
            fracs <- apply(cnts[, numeric.cols], 2, function(vec) vec / sum(vec))
            rownames(fracs) <- rownames(cnts)
            t(fracs)
          })
  colnames(frac.df)[1] <- "sample"

  frac.df[is.na(frac.df)] <- 0
  numeric.cols <- colnames(frac.df)[!(colnames(frac.df) %in% c("sample", "gene_biotype"))]
  frac.sums <- colSums(frac.df[, numeric.cols], na.rm=TRUE)
  frac.sums <- frac.sums[order(frac.sums, decreasing=TRUE)]
  print(frac.sums)
  max.biotypes.to.display <- 5
  other.biotypes <- names(frac.sums)[(max.biotypes.to.display+1):length(frac.sums)]
  fracs <- cbind(frac.df, other = rowSums(frac.df[, other.biotypes]))
  biotype.cols <- c(names(frac.sums)[1:max.biotypes.to.display],"other")
  #fracs <- as.data.frame(fracs)
  foo <- reshape2::melt(fracs[, c("sample", biotype.cols)])
  colnames(foo) <- c("sample", "biotype", "proportion")
  foo$proportion <- foo$proportion * 100
  foo$biotype <- factor(foo$biotype, biotype.cols)
  if(!is.null(sample.labels)) {
    sample.df <- data.frame(sample = names(sample.labels), label = as.character(sample.labels))
    foo <- merge(foo, sample.df)
    foo$sample <- factor(foo$sample)
  }
  g <- ggplot(data = foo, aes(x = biotype, y = proportion)) + geom_boxplot()
  if(is.null(sample.labels)) {
    g <- g + facet_wrap(~ sample, nrow=2)
  } else {
    g <- g + facet_wrap(label ~ sample)
  }
  g <- g + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  g <- g + xlab("Biotype") + ylab("Proportion")
  g
}

#' Create boxplots showing the contribution to total expression of the top n genes
#'  
#' @param mat An expression matrix
#' @param n.top The number of top genes (by frequency) to plot
#' @param highlight.genes A list of genes to highlight in box
#' @return A ggplot
plot.top.genes <- function(mat, n.top = 20, highlight.genes = NULL) {
  mat <- sweep(mat, 2, colSums(mat), "/")
  top.genes <- get.top.genes.matrix(mat, n.top = n.top)
  df <- melt(mat[top.genes,])
  colnames(df) <- c("gene", "spot", "value")
  df$gene <- factor(df$gene, levels = rev(top.genes))
  g <- ggplot() + geom_boxplot(data = df, aes(x = gene, y = 100 * value), fill = (scales::hue_pal())(n.top)[n.top:1]) + coord_flip()
  g <- g + ylab("% total count per cell") + xlab("")
  g <- g + theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))
  if(!is.null(highlight.genes)) {
    vec_fontface <- ifelse(levels(df$gene) %in% highlight.genes,"bold","plain")
    g <- g + theme(axis.text.y=element_text(face=vec_fontface))
  }
  g
}
