#' Wrapper function that calls all basic qc plotting funcitons
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' 
plot.qc.data<-function(dataset.data){
  # Plot genes and Umis for unfiltered data
  plot.basic.spatialplots(dataset.data,TRUE)
  # Plot genes and Umis for filtered data
  plot.basic.spatialplots(dataset.data,FALSE)
  # Plot biotypes per sample
  plot.biotypes.multiple.samples(dataset.data,'per_sample')
  # Plot biotypes per protocol (FFPE/FF)
  plot.biotypes.multiple.samples(dataset.data,'per_protocol')
  #Plot deconvolution results using unsupervised deconvolution STdeconvolution
  plot.deconvolution.results(dataset.data$cell_type_deconvolution$stdeconvolution.results,dataset.data,'per sample','ST')
  
  # In case there are single cell data and RCTD deconvolution has been calculated
  if ("correlation" %in% names(dataset.data[["cell_type_deconvolution"]])){
    # Plot RCTD results 
    plot.deconvolution.results(dataset.data$cell_type_deconvolution$rctd.results,dataset.data,'per sample','RCTD')
    # Plot correlation among RCTD and ST results
    plot.deconvolution.correlation.results(dataset.data)
  }
  
  # In case there are both FFPE and FF preserved samples
  if(length(unique(dataset.data$metadata$labels))==2){
    # Plot ffpe vs ff mean gene expression values
    plot.ffpe.ff.expression.density(dataset.data)
    # Plot GO enrichment for genes that exist only in FFPE samples
    plot.enirchment.results(dataset.data[["FFPE.only.enrichement.analysis"]][["GO.enrichment"]],dataset.data,'G0_FFPE')
    # Plot GO enrichment for genes that exist only in FFPE samples
    plot.enirchment.results(dataset.data[["FF.only.enrichement.analysis"]][["GO.enrichment"]],dataset.data,"GO_FF")
    # Plot KEGG pathways for genes that exist only in FFPE samples
    plot.enirchment.results(dataset.data[["FFPE.only.enrichement.analysis"]][["KEGG.analysis"]],dataset.data,'KEGG_FFPE')
    # Plot KEGG pathways for genes that exist only in FFPE samples
    plot.enirchment.results(dataset.data[["FF.only.enrichement.analysis"]][["KEGG.analysis"]],dataset.data,'KEGG_FF')
  }
}


#' Plot basic spatial plots including H&E staining, unfiltered gene/UMIs counts
#'
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @param filtered boolean variable indicating whether of not filtered values would be shown
plot.basic.spatialplots<-function(dataset.data, filtered=TRUE){
  datasets <- dataset.data$metadata$names
  names(datasets) <- dataset.data$metadata$names
  hne.plts <-   
    llply(datasets,
          .fun = function(dataset) {
            g<-plot.hne(dataset.data$filtered.objs[[dataset]], keep.invisible.legend = TRUE)
            g <- add.title.to.plot(g, paste0(dataset))
            png(paste0(dataset.data$plots_dir, "/", dataset, "-hne.png"))
            print(g)
            d <- dev.off()
            g
          })
  g.all <- plot_grid(plotlist = hne.plts)
  png(paste0(dataset.data$plots_dir, "/", "all-datasets-hne.png"))
  print(g.all)
  d <- dev.off()
  
  gene.plts <-   
    llply(datasets,
          .fun = function(dataset) {
            if (filtered==TRUE){
              g<-plot.features(dataset.data$filtered.objs[[dataset]], features=NULL, feature.names = NULL, include.hne = FALSE, include.umi.cnts = TRUE, include.feature.cnts = TRUE) 
              filter_prefix='filtrered'
            }else{
              g<-plot.features(dataset.data$unfiltered.objs[[dataset]], features=NULL, feature.names = NULL, include.hne = FALSE, include.umi.cnts = TRUE, include.feature.cnts = TRUE)
              filter_prefix='unfiltered'
            }
            g <- add.title.to.plot(g, paste0(dataset))
            png(paste0(dataset.data$plots_dir, "/", dataset,filter_prefix, "-umis-genes.png"),width=1080,height=680)
            print(g)
            d <- dev.off()
            g
          })
  g.all <- plot_grid(plotlist = gene.plts)
  if (filtered==TRUE){
    png(paste0(dataset.data$plots_dir, "/",'filtered', " all-datasets-umis-genes.png"),width=1480,height=480)
  }else{
    png(paste0(dataset.data$plots_dir, "/","ufiltered", " all-datasets-umis-genes.png"),width=1480,height=480)
  }
  print(g.all)
  d <- dev.off()
  
  temp_metadata<-c('num.umis','num.nonzero.genes','saturation')
  temp_metadata_ylab_names<-c("Number of UMIs","Number of Genes\n(with counts > 0)","Saturation","Fraction of mapped reads")
  if (filtered==TRUE){
    filter_prefix='filtrered'
    temp_filtered=dataset.data$filtered.metadata[dataset.data[["filtered.metadata"]][["spot_type"]]=='tissue',]
  }else{
    filter_prefix='unfiltrered'
    temp_filtered=dataset.data$filtered.metadata
  }
  wrap_var<-'per_protocol'
  plots<-
    llply(1:length(temp_metadata),
          .fun=function(i){
            g<-plot.metadata.feature(temp_filtered,temp_metadata[i],temp_metadata_ylab_names[i],wrap_var)
            g<- g + theme(legend.position= "bottom" )
            png(paste0(dataset.data$plots_dir, "/",wrap_var," ", temp_metadata[i]," ",filter_prefix,".png"), width = 2 * 480)
            print(g)
            d <- dev.off()
            g
          })
  g<-plot.metadata.features(plots)
  png(paste0(dataset.data$plots_dir, "/", filter_prefix," all-num-umis-and-genes.png"), width = 2 * 480)
  print(g)
  d <- dev.off()
}

#' Plot deconvolution correlation results
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
plot.deconvolution.correlation.results<-function(dataset.data){
  datasets <- dataset.data$metadata$names
  names(datasets) <- dataset.data$metadata$names
  hne.plts <-   
    llply(datasets,
          .fun = function(dataset) {
            png(paste0(dataset.data$plots_dir, "/", dataset," Correlation_matrix_RCTD_vs_STDeconvole.png"), width = 2 * 480)
            corrplot(dataset.data[["cell_type_deconvolution"]][["correlation"]][[dataset]],title=paste0(dataset),mar = c(0,0,1,0))
            dev.off()
          })
}

#' Plot deconvolution results
#' 
#' @param deconvolution.df A dataframe where each row is a spot for each of the samples that we have deconvoluted their values
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @param wrap_var  A variable specifying whether our plots would be seperated through protocol or not
#' 
plot.deconvolution.results<-function(deconvolution.df,dataset.data,wrap_var='per sample',prefix){
  temp_metadata<-colnames(deconvolution.df)[!colnames(deconvolution.df) %in% c("label","orig.ident","spot_type")]
  temp_metadata_ylab_names<-c("Prc of cell type")
  plots<-
    llply(1:length(temp_metadata),
          .fun=function(i){
            g<-plot.metadata.feature(deconvolution.df,temp_metadata[i],temp_metadata_ylab_names,wrap_var)
            g<- g + theme(legend.position= "bottom" )
            png(paste0(dataset.data$plots_dir, "/",temp_metadata[i]," prc.png"), width = 2 * 480,height=480)
            print(g)
            d <- dev.off()
            g
          })
  
  g <- plot.population.fractions.across.samples(deconvolution.df, id.cols = c("orig.ident", "label","spot_type"), sample.col = c("orig.ident"), fill = "label", linetype = "spot_type") + labs(fill="Protocol") + theme(text = element_text(size = 20))
  png(paste0(dataset.data$plots_dir, "/",prefix," all cell types prc.png"), width = 2 * 480)
  print(g)
  d <- dev.off()
  # g<-plot.metadata.features(plots)
  # png(paste0(dataset.data$plots_dir, "/all cell types prc.png"), width = 2 * 480,height=480)
  # print(g)
  # d <- dev.off()
}

#' Create a plot grid out of multiple plots
#' 
#' @param plots Grid of plots
#' 
plot.metadata.features<-function(plots){
  g.legend <- get_legend(plots[[1]])
  for(i in 1:length(plots)){
    plots[[i]] <- plots[[i]] + theme(legend.position = "none")
  }
  g <- plot_grid(plotlist=plots, ncol = length(plots), labels = "AUTO")
  g <- plot_grid(g, g.legend, nrow = length(plots), rel_heights=c(5,1))
  g
}

#' Create a plot of a specific metadata feature
#' 
#' @param all.unfiltered.meta A dataframe where each row is a spot for each of the samples and columns are the calculated metadata
#' @param temp_feature The specific metadata that we want to plot
#' @param wrap_var  A variable specifying whether our plots would be seperated through protocol or not
#' 
plot.metadata.feature<-function(all.unfiltered.meta,temp_feature,temp_ylab, wrap_var=NULL){
  g <- ggplot(data = all.unfiltered.meta, aes(x = orig.ident, y = get(temp_feature), fill = label, linetype = spot_type))
  if(identical(wrap_var,'per_protocol')) {
    g <- g + facet_wrap(~ label, scales="free")
  }
  g <- g + ggtitle(temp_feature)
  g <- g + theme(plot.title = element_text(hjust = 0.5,size=16, face = "bold"))
  g <- g + geom_violin() + scale_y_continuous(trans="log2")
  g <- g + labs(fill="Protocol") + labs(linetype="Spot") + ylab(temp_ylab) + xlab("")
  g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  g <- g + theme(legend.position = "top", text = element_text(size = 15), axis.text.y = element_text(size=16, face = "bold"))
  g
}

#' Plot biotypes
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @param wrap_var  A variable specifying whether our plots would be seperated through protocol or not
#' @param title Character variable used as title for the ggplot object
plot.biotypes.multiple.samples<-function(dataset.data,wrap_var, title=NULL){
  frac.df<-dataset.data$biotypes
  frac.df[is.na(frac.df)] <- 0
  numeric.cols <- colnames(frac.df)[!(colnames(frac.df) %in% c("sample", "gene_biotype"))]
  frac.sums <- colSums(frac.df[, numeric.cols], na.rm=TRUE)
  frac.sums <- frac.sums[order(frac.sums, decreasing=TRUE)]
  max.biotypes.to.display <- 5
  other.biotypes=NULL
  if (length(frac.sums)>max.biotypes.to.display){
    other.biotypes <- names(frac.sums)[(max.biotypes.to.display+1):length(frac.sums)]
    if (length(other.biotypes)>1){
      fracs <- cbind(frac.df, other = rowSums(frac.df[, other.biotypes]))
    } else {
      fracs <- cbind(frac.df, other = frac.df[, other.biotypes])
    }
    biotype.cols <- c(names(frac.sums)[1:max.biotypes.to.display],"other")
  }else{
    fracs <- cbind(frac.df, other = frac.df[, other.biotypes])
    biotype.cols <- c(names(frac.sums))
  }
  #fracs <- as.data.frame(fracs)
  foo <- reshape2::melt(fracs[, c("sample", biotype.cols)])
  colnames(foo) <- c("sample", "biotype", "proportion")
  foo$proportion <- foo$proportion * 100
  foo$biotype <- factor(foo$biotype, biotype.cols)
  if(!is.null(dataset.data$metadata$labels)) {
    sample.df <- data.frame(sample = dataset.data$metadata$names, label = as.character(dataset.data$metadata$labels))
    foo <- merge(foo, sample.df)
    foo$sample <- factor(foo$sample)
  }
  g <- ggplot(data = foo, aes(x = biotype, y = proportion)) + geom_boxplot()
  if(is.null(dataset.data$metadata$labels)) {
    g <- g + facet_wrap(~ sample, nrow=2)
  } else {
    if(wrap_var=='per_protocol'){
      g <- g + facet_wrap(~ label)
    }else if (wrap_var=='per_sample'){
      g <- g + facet_wrap(~ label ~ sample)
    }
  }
  g <- g + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  g <- g + xlab("Biotype") + ylab("Proportion")
  g<- g + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
  png(paste0(dataset.data$plots_dir, "/ ",wrap_var,"-biotypes.png"), width = 2 * 480)
  print(g)
  d <- dev.off()
}

#' Plot enirchment results
#' 
#' @param Kegg_enrichment Enrichment result from the clusterprofiler package
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @param analysis_set Character variable used as title for the ggplot object an file description
#' 
plot.enirchment.results<-function(enrichment,dataset.data,file_prefix){
  g<-  dotplot(enrichment,showCategory = 15)
  png(paste0(dataset.data$plots_dir, "/", "dotplot_","_" ,file_prefix ,".png"), width = 800, height = 800)
  print(g)
  d <- dev.off()
  p.go <- pairwise_termsim(enrichment)
  g.go <- emapplot(p.go)
  png(paste0(dataset.data$plots_dir, "/", "emmaplot_", file_prefix ,".png"), width = 1000, height = 1000)
  print(g.go + theme(text = element_text(size = 20)))
  d <- dev.off()
}

#' Plot cell type deconvolution

#' Plot density plot FFPE vs FF genes mean expression
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
plot.ffpe.ff.expression.density<-function(dataset.data){
  m <- merge(dataset.data[["gene.exrpession"]][["gene.expr.ff"]], dataset.data[["gene.exrpession"]][["gene.expr.ffpe"]], all = TRUE)
  # See https://github.com/daattali/ggExtra/issues/27
  g <- ggplot(data = m, aes(x = log2(1+mean.expr.frozen), y = log2(1+mean.expr.ffpe)))
  g <- g + xlab("Frozen Log2(1 + Total UMIs Across Spots)") + ylab("FFPE Log2(1 + Total UMIs Across Spots)")
  g <- g + geom_point(col="transparent")+geom_hex()+labs(fill="UMIs")
  g <- g + theme(text = element_text(size = 15))
  g.leg <- get_legend(g)
  g <- ggMarginal(g + theme(legend.position = "none"), type="histogram")
  g <- plot_grid(g, g.leg, rel_widths=c(10,1))
  
  png(paste0(dataset.data$plots_dir, "/", "ffpe-vs-frozen-scatter.png"), width = 1 * 480)
  print(g)
  d <- dev.off()
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
plot.features <- function(obj, features, feature.names = NULL, include.hne = FALSE, include.umi.cnts = FALSE, include.feature.cnts = FALSE, ...) {
  plts <- NULL
  if(!is.null(features) && (length(features) > 0)) {
    if(is.null(feature.names)) { feature.names <- features }
    indices <- 1:length(features)
    names(indices) <- features
    plts <- lapply(indices, function(i) plot.spatial(obj, features = c(features[i]), legend.name = feature.names[i], rescale.legend = FALSE))
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
  plot_grid(plotlist = plts, ...)
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
  g <- g + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15), legend.position = "right")
  # See https://stackoverflow.com/questions/45998396/unset-existing-scale-fill-discrete-in-ggplot2-or-suppress-message-for-new-scale
  # for the following code, which suppresses 
  # "Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale."
  # g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0fk', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  i <- which(sapply(g$scales$scales, function(x) 'fill' %in% x$aesthetics))
  g$scales$scales[[i]] <- NULL
  if(rescale.legend) {
    g <- g + scale_fill_gradientn(name = legend.name, labels = function(x) { sprintf('%.0f', x/1000) }, colours = Seurat:::SpatialColors(n = 100))
  } else {
    g <- g + scale_fill_gradientn(name = legend.name, colours = Seurat:::SpatialColors(n = 100))
  }
  # g <- g + theme(legend.text = element_text(angle = 45, vjust = 1, hjust=1))
  g <- g + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  g.leg <- get_legend(g) # + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  plot_grid(g + theme(legend.position = "none"), g.leg, nrow = 1, rel_widths = c(7,3))
  # g
}