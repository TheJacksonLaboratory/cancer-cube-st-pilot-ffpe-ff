#' Wrapper function that calls all basic qc functions
#' 
#' @param filename A configuration file including all the necessary information for loading the data of that dataset 
#' @param output_dir A base directory where the config .csv file will be saved and the analysis files 
#' @return dataset.data A list with all the important metadata and processed entities after initial QC
#' 
#' 
calculate.qc.functions<-function(filename){
  # Define output directory as the folder where the config file is located 
  output_dir <- paste0(dirname(filename),'/')
  # Load the metadata of the dataset and the filtered and unfiltered data of each of the samples of the dataset
  dataset.data <- load.init.file(filename, output_dir)
  # Calculate for each spot of each sample of the dataset useful metadata (described in detail in the retrieve.metrics function below)
  dataset.data$filtered.metadata <- retrieve.metrics(dataset.data)
  # Calculate for each spot of each sample of the dataset the biotypes that the genes that expressed are asosciated with (based on the species that the dataset is derived from) 
  dataset.data$biotypes <- calculate.biotypes.per.spot(dataset.data)
  # Calculate the expression value of each of the genes that is expressed in the dataset, the mean expression of the genes in ffpe/ff samples
  dataset.data$gene.exrpession <- gene.exrpession(dataset.data)
  
  # Normalize the data based on the following normalization strategies
  normalization_strategy <- c("SCT","SimpleNorm","cpm","logcpm")
  for(i in normalization_strategy){
    dataset.data$filtered.objs <- dataset.normalization(dataset.data,i)
  }
  
  # Calculate cell type deconvolution
  if(!is.na(dataset.data$metadata$single_cell_expr[1]) & !is.na(dataset.data$metadata$single_cell_cell_types[1])){
    # In case there are single.cell data calculate in a supervised the cell types percentages of each spot, based on RCTD deconvolution
    dataset.data$cell_type_deconvolution$rctd.results <- rctd.deconvolution(dataset.data)
    # Define the number of expected cell types ( input argument for STdeconvolution) - here defined as the same as the RCTD in order to calculate their correlation
    opt_val=length(dataset.data[["cell_type_deconvolution"]][["rctd.results"]])-3
  }else{
    # Define the number of expected cell types ( input argument for STdeconvolution) - here defined as the minimun of the calculated LDA models
    opt_val='min'
  }
  
  # Calculate ST cell type deconvolution, in an unsupervised way regardless the existance of single-cell data
  dataset.data$cell_type_deconvolution$stdeconvolution.results <- stdeconvolve.deconvolution(dataset.data,opt_val)
  if(!is.na(dataset.data$metadata$single_cell_expr[1]) & !is.na(dataset.data$metadata$single_cell_cell_types[1])){
    dataset.data$cell_type_deconvolution$correlation <- deconvolution.correlation(dataset.data)
  }
  
  # In case there are samples derived from FFPE and FF preseration protocol
  if(length(unique(dataset.data$metadata$labels))==2){
    #Calculate GO enrichment and KEGG pathways for FFPE only genes
    dataset.data$FFPE.only.enrichement.analysis <- GO.KEGG.enirhcment.analysis(dataset.data$gene.exrpession$FFPE.only.genes,dataset.data$metadata$species[1])
    #Calculate GO enrichment and KEGG pathways for FF only genes
    dataset.data$FF.only.enrichement.analysis <- GO.KEGG.enirhcment.analysis(dataset.data$gene.exrpession$FF.only.genes,dataset.data$metadata$species[1])
  }
  dataset.data
}


#' Load a .csv file for all the necessary information needed for the pre-proceesing
#' of a specific dataset
#' 
#' @param filename A file including all the necessary information directory where the .csv file will be saved#' 
#' @param output_dir A base directory where the config .csv file will be saved and the analysis files.
#'                   If null (default), then output analysis and plots dirs are not created and neither is .csv config file 
#' @return temp_list A list with all the important metadata and filtered and unfiltered values as they 
#' are described here: https://satijalab.org/seurat/reference/load10x_spatial
#' 
load.init.file <- function(filename, output_dir = NULL) {
  metadata <- read.csv(file = filename, stringsAsFactors = FALSE)
  metadata$tissue_type <-as.factor(basename(dirname(filename)))
  spaceranger_dirs <- metadata$dataset_names
  names(spaceranger_dirs) <- metadata$names
  filtered.objs <- create.visium.seurat.objects(spaceranger_dirs, filter.spots = TRUE)
  unfiltered.objs <- create.visium.seurat.objects(spaceranger_dirs, filter.spots = FALSE)
  temp_list<-list(metadata=metadata,filtered.objs=filtered.objs,unfiltered.objs=unfiltered.objs)

  if(!is.null(output_dir)) {
    temp_list$base_dir<-output_dir
    temp_list$analysis_dir <- paste0(output_dir, '/analysis/')
    dir.create(temp_list$analysis_dir, showWarnings = FALSE, recursive = TRUE)
    temp_list$plots_dir <- paste0(output_dir, '/plots/')
    dir.create(temp_list$plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Write csv file for automatic ST pipeline 
    write.automated.csv.file(output_dir,temp_list$metadata$names, temp_list$metadata$dataset_names, temp_list$metadata$species)
  }
  temp_list
}

#'Create a dataframe with all basic metric statistics as they are retrieved from 
#'spaceranger metadata files and return the same object with a field metadata_retrieved 
#'containing a dataframe with rows each spot and columns:
#'  E: The number of exonic reads in the spot
#'  I: The number of intronic reads in the spot
#'  N: The number of intergenic reads in the spot
#'  tot.mapped.reads: The number of mapped reads in the spot
#'  tot.unmapped.reads: The number of unmapped reads in the spot
#'  tot.reads: tot.mapped.reads + tot.unmapped.reas
#'  E.frac.mapped: E / tot.mapped.reads
#'  I.frac.mapped: I / tot.mapped.reads
#'  N.frac.mapped: N / tot.mapped.reads
#'  num.reads: total number of reads to barcode
#'  num.umis: number of unique reads (UMIs) to barcode
#'  num.nonzero.genes: number of genes to barcode with one or more reads
#'  saturation: (num.reads - num.unis) / num.reads
#'
#'
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return all.filtered.meta A dataframe with all meta data 
#'
retrieve.metrics<-function(dataset.data){
  nms <- dataset.data$metadata$dataset_names
  names(nms) <- dataset.data$metadata$names
  
  align.metrics <-
  llply(nms, .parallel = FALSE,
        .fun = function(nm) {
          alignment.metric.file <- paste0(dataset.data$analysis_dir, "/", names(nm), "-alignment-metrics.csv")
          if(!file.exists(alignment.metric.file)) {
            bam.file <- paste0(nm, "/possorted_genome_bam.bam")
            df <- get.per.spot.alignment.metrics(bam.file)
            write.table(file=alignment.metric.file, df, row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)
          }
          as.data.frame(fread(alignment.metric.file))
        })
  
  filtered.saturation.metrics <-
    llply(nms, .parallel = TRUE,
          .fun = function(nm) {
            mol.info.file <- paste0(nm, "molecule_info.h5")
            get.per.spot.saturation.stats(mol.info.file, filter = TRUE)
          })
  
  print('finished saturation')
  temp_list_<-list()
  for(nm in names(nms)){
    # align.meta<-list()
    align.meta <- align.metrics[[nm]]
    # Spot barcode is in CB column. Make that the row name and remove CB.
    align.meta <- align.meta[!is.na(align.meta$CB),]
    rownames(align.meta) <- align.meta$CB
    align.meta <- align.meta[, !(colnames(align.meta) == "CB")]
    
    seq.meta <- filtered.saturation.metrics[[nm]]
    # Add "-1" postfix to barcode
    rownames(seq.meta) <- paste0(rownames(seq.meta), "-1")
    m <- merge(align.meta, seq.meta, by = "row.names", all = TRUE)
    m[is.na(m)] <- 0
    rownames(m) <- m$Row.names
    m <- m[, !(colnames(m) == "Row.names")]
    m$map.frac <- m$tot.mapped.reads / m$tot.reads
    temp<-intersect(dataset.data[["unfiltered.objs"]][[nm]]@assays[["Spatial"]]@data@Dimnames[[2]],rownames(m))
    m=m[temp,]
    m$spot_type<-dataset.data[["unfiltered.objs"]][[nm]]@meta.data[["spot_type"]][match(temp,dataset.data[["unfiltered.objs"]][[nm]]@assays[["Spatial"]]@data@Dimnames[[2]])]
    temp_list_[[nm]]<-m
  }
  all.filtered.meta <- ldply(temp_list_, .fun = function(obj) obj)
  colnames(all.filtered.meta)[colnames(all.filtered.meta) == ".id"] ="orig.ident"
  sample.df <- data.frame(sample = dataset.data$metadata$names, label = as.character(dataset.data$metadata$labels))
  all.filtered.meta <- merge(all.filtered.meta, sample.df, by.x = c("orig.ident"), by.y = c("sample"))
  all.filtered.meta$orig.ident <- factor(all.filtered.meta$orig.ident)
  all.filtered.meta
}


#' Calculate the different biotypes the genes of a dataset are expressing having as
#' input an object with counts and metadata and return
#' a dataframe with the biotypes of each spot of each sample.
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return biotypes A dataframe where each row is a spot and columns the biotypes of the genes expressed in these spots
#'
calculate.biotypes.per.spot<-function(dataset.data){
  if (dataset.data[["metadata"]][["species"]][1] == 'Human'){
    gene_db = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  } else {
    gene_db = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  biotypes <-
    ldply(dataset.data[["filtered.objs"]], .parallel = FALSE,
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
  colnames(biotypes)[1] <- "sample"
  biotypes
}

#' Calculation of mean expression value of each gene either in tissue or background
#' 
#'
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return gene.expression A list containing info for the expression of each gene, mean ffpe/ff epxression, if applicable
#' 
gene.exrpession<-function(dataset.data, filtering_tissue_spots=TRUE){
  datasets <- dataset.data$metadata$names
  names(datasets) <- dataset.data$metadata$names
  all.genes <-
    ldply(datasets,
          .fun = function(dataset) {
            obj <- dataset.data$unfiltered.objs[[dataset]]
            mat <- as.matrix(GetAssayData(obj, slot="counts", assay="Spatial"))
            if (filtering_tissue_spots==TRUE){
              mat <- mat[,dataset.data[["unfiltered.objs"]][[dataset]]@meta.data[["tissue"]]==1]
            }
            temp_val=as.vector(rowSums(mat))
            data.frame(gene = rownames(mat), expr.sum = temp_val/sum(temp_val) )
          })
  colnames(all.genes)[colnames(all.genes) == ".id"] ="orig.ident"
  sample.df <- data.frame(sample = dataset.data$metadata$names, label = as.character(dataset.data$metadata$labels), spot_type='tissue')
  all.genes <- merge(all.genes, sample.df, by.x = c("orig.ident"), by.y = c("sample"))
  
  gene.expr.ffpe <-
    ddply(subset(all.genes, orig.ident %in% names(datasets)[dataset.data$metadata$labels=="FFPE"]),
          .variables = c("gene"),
          .fun = function(df) {
            data.frame(mean.expr.ffpe = mean(df$expr.sum))
          })
  gene.expr.frozen <-
    ddply(subset(all.genes, orig.ident %in% names(datasets)[dataset.data$metadata$labels=="Frozen"]),
          .variables = c("gene"),
          .fun = function(df) {
            data.frame(mean.expr.frozen = mean(df$expr.sum))
          })
  ffpe.ff.corr<-0
  FFPE.only.genes<-0
  FF.only.genes<-0
  if (nrow(gene.expr.frozen)>0 & nrow(gene.expr.ffpe)>0){
    m <- merge(gene.expr.frozen, gene.expr.ffpe, all = TRUE)
    m[is.na(m)] <- 0
    cor_m<-m
    # For correlation we keep only targeted genes in FFPE samples
    cor_m<-cor_m[cor_m$gene %in% dataset.data[["filtered.objs"]][[dataset.data$metadata$names[dataset.data$metadata$labels=='FFPE'][1]]]@assays[["Spatial"]]@data@Dimnames[[1]],]
    ffpe.ff.corr<-cor(cor_m$mean.expr.frozen,cor_m$mean.expr.ffpe)
    FFPE.only.genes<-m$gene[m$mean.expr.ffpe>0 & m$mean.expr.frozen==0]
    FF.only.genes<-m$gene[m$mean.expr.ffpe==0 & m$mean.expr.frozen>0]
  }
  gene.expressions<-list(all.genes=all.genes,gene.expr.ffpe=gene.expr.ffpe,gene.expr.ff=gene.expr.frozen, correlation=ffpe.ff.corr, FF.only.genes=FF.only.genes,FFPE.only.genes=FFPE.only.genes)
}

#' Deconvolution correlation
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
deconvolution.correlation<-function(dataset.data){
  datasets <- dataset.data$metadata$names
  names(datasets) <- datasets
  main.rctds <- 
    llply(datasets, .parallel = FALSE,
          .fun = function(dataset) {
            obj <- dataset.data$filtered.objs[[dataset]]
            cd <- GetAssayData(obj, assay="Spatial", slot="counts")
            counts <- cleanCounts(counts = cd,
                                  min.lib.size = 100,
                                  min.reads = 1,
                                  min.detected = 1,
                                  verbose = TRUE)
            corpus <- restrictCorpus(counts, removeAbove=0.98, removeBelow = 0.05, nTopOD = 1000)
            rds.file.STDe <- paste0(dataset.data$analysis_dir, "/",dataset, "-STdeconvolve.rds")
            if(!file.exists(rds.file.STDe)) {
              print('First perform ST deconvolution')
            } else  {
              cat(paste0("Reading ", rds.file.STDe, "\n"))
              optLDA <- readRDS(rds.file.STDe)
            }
            myRCTD <-getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
            myRCTD
            
            rds.file <- paste0(dataset.data$analysis_dir, "/", dataset, "-rctd-main.rds")
            if(!file.exists(rds.file)) {
              print('First perform RCTD deconvolution')
            } else {
            cat(paste0("Reading ", rds.file, "\n"))
            myRCTD2 <- readRDS(rds.file)
            }
            gc()
            temp <- intersect(rownames(myRCTD$theta),myRCTD2@results[["weights"]]@Dimnames[[1]])
            wqer <- cor(as.matrix(myRCTD2@results[["weights"]][temp,]),as.matrix(myRCTD[["theta"]][temp,]))
            wqer
          })
  main.rctds
}

#' Calculation for each spot the percentage of cells that belong to a specific cell type
#' using RCTD method (doi:10.1038/s41587-021-00830-w)
#'
#' @param counts_sc A dataframe where columns are the cell barcodes and rows the raw counts of the expressed genes 
#' @param cell_types  large factor where the cell type of each cell is described, with names the barcodes as they are defined in the param above
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return main.stdeconvolutions A list with each element of the list shows the results of each sample of the dataset
#' 
stdeconvolve.deconvolution<-function(dataset.data,opt_val='min'){
  datasets <- dataset.data$metadata$names
  names(datasets) <- datasets
  main.rctds <- 
    llply(datasets, .parallel = FALSE,
          .fun = function(dataset) {
            obj <- dataset.data$filtered.objs[[dataset]]
            cd <- GetAssayData(obj, assay="Spatial", slot="counts")
            counts <- cleanCounts(counts = cd,
                                  min.lib.size = 100,
                                  min.reads = 1,
                                  min.detected = 1,
                                  verbose = TRUE)
            corpus <- restrictCorpus(counts, removeAbove=0.98, removeBelow = 0.05, nTopOD = 1000)
            rds.file.STDe <- paste0(dataset.data$analysis_dir, "/",dataset, "-STdeconvolve.rds")
            if(!file.exists(rds.file.STDe)) {
              if(opt_val=='min'){
                ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(5, 11, by = 2),perc.rare.thresh = 0.05, plot=TRUE, verbose=TRUE)
                optLDA <- optimalModel(models = ldas, opt = opt_val)
              } else{
                ldas <- fitLDA(t(as.matrix(corpus)), Ks = opt_val,perc.rare.thresh = 0.05, plot=TRUE, verbose=TRUE)
                optLDA <- optimalModel(models = ldas, opt = opt_val)
              }
              saveRDS(optLDA, rds.file.STDe)
            } else  {
              cat(paste0("Reading ", rds.file.STDe, "\n"))
              optLDA <- readRDS(rds.file.STDe)
            }
            myRCTD <-getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
            myRCTD
          })
  all.deconvolution <- ldply(main.rctds, .fun = function(obj) as.data.frame(obj[["theta"]]))
  colnames(all.deconvolution)[colnames(all.deconvolution) == ".id"] ="orig.ident"
  sample.df <- data.frame(sample = dataset.data$metadata$names, label = as.character(dataset.data$metadata$labels), spot_type='tissue')
  all.deconvolution <- merge(all.deconvolution, sample.df, by.x = c("orig.ident"), by.y = c("sample"))
  all.deconvolution$orig.ident <- factor(all.deconvolution$orig.ident)
  all.deconvolution
}

#' Calculation for each spot the percentage of cells that belong to a specific cell type
#' using RCTD method (doi:10.1038/s41587-021-00830-w)
#'
#' @param counts_sc A dataframe where columns are the cell barcodes and rows the raw counts of the expressed genes 
#' @param cell_types  large factor where the cell type of each cell is described, with names the barcodes as they are defined in the param above
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return main.rctds A list with each element of the list shows the results of each sample of the dataset
#' 
rctd.deconvolution<-function(dataset.data){
  
  datasets <- dataset.data$metadata$names
  names(datasets) <- datasets
  main.rctds <- 
    llply(datasets, .parallel = FALSE,
          .fun = function(dataset) {
            print(dataset)
            rds.file <- paste0(dataset.data$analysis_dir, "/", dataset, "-rctd-main.rds")
            if(!file.exists(rds.file)) {
              obj <- dataset.data$filtered.objs[[dataset]]
              st.counts <- GetAssayData(obj, assay="Spatial", slot="counts")
              st.coords <- obj[[]][, c("col", "row")]
              colnames(st.coords) <- c("x","y")
              nUMI <- colSums(st.counts) # In this case, total counts per pixel is nUMI
              puck <- SpatialRNA(st.coords, st.counts, nUMI)
              
              # Create a reference for _each_ sample. Why? Because the Reference is
              # dependent on a single-cell RNA-seq reference (sc_mat) and we want
              # to ensure that those genes are in the spots. If they weren't, we
              # may have genes selected as markers that aren't actually in the ST
              # data.
              
              # retrieve the targeted genes, asa they are mentioned in the fitlered object
              gt<-dataset.data[["filtered.objs"]][[dataset.data$metadata$names[dataset.data$metadata$labels=='FFPE'][1]]]@assays[["Spatial"]]@data@Dimnames[[1]]
              
              # Read single-cell expression file
              counts_sc_full <- fread(dataset.data$metadata$single_cell_expr[match(dataset,dataset.data$metadata$names)])
              foo <- as.data.frame(counts_sc_full)
              rownames(foo) <- foo$index
              counts_sc <- t(foo[,-1])
              # Read single-cell cell type annotation file
              cell_types <- read.csv(dataset.data$metadata$single_cell_cell_types[match(dataset,dataset.data$metadata$names)])
              cell_types_main<- setNames(cell_types[[2]], cell_types[[1]])
              cell_types_main <- as.factor(cell_types_main) # convert to factor data type
              
              sc_mat <- counts_sc[rownames(counts_sc) %in% gt,]
              cell_types <- cell_types_main[colnames(sc_mat)]
              nUMI <- colSums(sc_mat)
              reference <- Reference(sc_mat, cell_types, nUMI, n_max_cells = ncol(sc_mat) + 1)
              
              max_cores <- min(10, num.cores - 1)
              # Note that we have keep_reference = TRUE here. Without it, the default is
              # keep_reference = FALSE, which would ignore the Reference we have 
              # intentionally created above.
              myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, keep_reference = TRUE)
              myRCTD <- suppressPackageStartupMessages(run.RCTD(myRCTD, doublet_mode = 'full'))
              saveRDS(myRCTD, rds.file)
            }
            cat(paste0("Reading ", rds.file, "\n"))
            myRCTD <- readRDS(rds.file)
            gc()
            myRCTD
          })
  all.deconvolution <- ldply(main.rctds, .fun = function(obj) as.data.frame(obj@results[["weights"]]))
  colnames(all.deconvolution)[colnames(all.deconvolution) == ".id"] ="orig.ident"
  sample.df <- data.frame(sample = dataset.data$metadata$names, label = as.character(dataset.data$metadata$labels), spot_type='tissue')
  all.deconvolution <- merge(all.deconvolution, sample.df, by.x = c("orig.ident"), by.y = c("sample"))
  all.deconvolution$orig.ident <- factor(all.deconvolution$orig.ident)
  all.deconvolution
}

#' Calculate GO and Gene Ontology and KEGG anaylsis
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @return enrichment.analysis A list containing the GO, KEGG results
GO.KEGG.enirhcment.analysis<-function(geneset,species){
  if (species == 'Human'){
    orgdb_temp = "org.Hs.eg.db"
    temp_organism='hsa'
  } else {
    orgdb_temp = "org.Mm.eg.db"
    temp_organism='mmu'
  }
  GO.enrichment <- enrichGO(geneset, OrgDb=orgdb_temp, keyType= 'SYMBOL', ont = "ALL", pvalueCutoff=0.1, qvalueCutoff = 0.1)
  
  gene.df <- bitr(geneset, fromType = "SYMBOL", toType = c("ENTREZID" ), OrgDb = orgdb_temp)
  KEGG.analysis <- enrichKEGG(gene = gene.df$ENTREZID , organism = temp_organism, pvalueCutoff = 0.1, qvalueCutoff = 0.1)
  
  enrichment.analysis<-list(GO.enrichment=GO.enrichment, KEGG.analysis=KEGG.analysis)
}

#' Dataset normialization
#' 
#' @param dataset.data A list containing Seurat objects and metadata information for the dataset
#' @param normalization_strategy A character vector indicating which normalization strategies should be applied to the dataset
dataset.normalization<-function(dataset.data, norm_assay='SCT'){
  if(norm_assay=='SCT'){
    # Perform normalization of counts using SCTransform (which fits a negative binomial to the count data)
    # This adds a new assay "SCT" to the returned Seurat object. That assay has slots "counts" (the corrected counts),
    # "data" (log1p(counts)), and "scale.data" (the "pearson residuals").
    # These would be accessible via GetAssayData(filtered.objs[[1]], slot = "data", assay = "SCT")
    # Running with default parameters (the most interesting of which look to be variable.features.n and/or variable.features.rv.th)
    dataset.data$filtered.objs <-
      llply(dataset.data$filtered.objs,
            .fun = function(obj) {
              obj <- subset(obj, nCount_Spatial>0)
              suppressWarnings(SCTransform(obj, assay = "Spatial", verbose = FALSE))
            })
  }else if(norm_assay=='logcpm'){
    # Add normalized (log) Counts Per Million (CPM) data to the filtered objects
    dataset.data$filtered.objs <-
      llply(dataset.data$filtered.objs,
            .fun = function(obj) {
              cpm.mat <- cpm(GetAssayData(obj, slot="counts", assay="Spatial"), log = TRUE)
              obj[["logcpm"]] <- CreateAssayObject(counts = cpm.mat)
              obj
            })
  }else if(norm_assay=='cpm'){
    # Add normalized (non-log) Counts Per Million (CPM) data to the filtered objects
    dataset.data$filtered.objs <-
      llply(dataset.data$filtered.objs,
            .fun = function(obj) {
              cpm.mat <- cpm(GetAssayData(obj, slot="counts", assay="Spatial"), log = FALSE)
              obj[["cpm"]] <- CreateAssayObject(counts = cpm.mat)
              obj
            })
  }else if(norm_assay=='SimpleNorm'){
    # Add normalized , with 'NormalizeData' function to the filtered objects
    dataset.data$filtered.objs <-
      llply(dataset.data$filtered.objs,
            .fun = function(obj) {
              cpm <- NormalizeData(obj, assay = "Spatial", verbose = FALSE,normalization.method='RC', scale.factor=1e6)
              obj[["SimpleNorm"]] <- CreateAssayObject(counts = as.matrix(cpm@assays[["Spatial"]]@data))
              obj
            })
  }
}

#' Create a .csv file for the automatic analysis of ST data, as it is described 
#' here: https://github.com/sdomanskyi/spatialtranscriptomics
#' 
#' @param base_dir A base directory where the .csv file will be saved
#' @param species_sample Vector describing the species of each sample
#' 
#' 
#' # Write csv file for automatic ST pipeline 
write.automated.csv.file <- function(base_dir,names_, dataset_names, species_sample) {
  temp_var2<- replicate(length(dataset_names), "")
  ST_automated<- data.frame (sample_id= names_, species= species_sample, st_data_dir=dataset_names, sc_data_dir= temp_var2 )
  dir.create(file.path(base_dir, '/automatic_analysis/'), recursive=TRUE)
  temp <- temp <- strsplit(base_dir,"/")
  curr_datatypes <- tail(temp[[1]], n=1)
  write.csv(ST_automated,paste0(base_dir,'/automatic_analysis/',curr_datatypes,'.csv'), row.names = FALSE ,quote=FALSE)
}

#' Create a Seurat object for a Visium 10X dataset.
#'  
#' @param sample.name A string name for the sample / dataset.
#' @param spaceranger.dir A string directory holding the files filtered_feature_bc_matrix.h5 and
#'                        raw_feature_bc_matrix.h5 and the subdirectory spatial. This is likely
#'                        the outs directory created by spaceranger count or the equivalent.
#' @param filter.spots A boolean indicating whether to include all spots (filter.spots = FALSE) or only
#'                     those spots overlapping the tissue (filter.spots = TRUE)
#' @return A Seurat object, annotated with the position and status (spot_type = "tissue" or "background") for each spot.
create.visium.seurat.object <- function(sample.name, spaceranger.dir, filter.spots = TRUE) {
  # load filtered_feature_bc_matrix.h5 with filter.matrix = TRUE or
  #      raw_feature_bc_matrix.h5 with filter.matrix = FALSE
  # filtered_feature_bc_matrix.h5 should only contain data from spots overlaying tissue,
  # though I can't find this documented conclusively.
  # filter.matrix = TRUE definitely filters coordinates in the Seurat object
  # according to those that overlap the tissue (i.e., have tissue == 1 in the
  # tissue_positions_list.csv file -- see code for Read10X_Image, which is called
  # by Load10X_Spatial) and for FFPE samples, non-targeted genes (https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices)
  
  # See this issue for how to set up these individual objects such that we can merge them:
  # https://github.com/satijalab/seurat/issues/3732
  # Namely, set the slice and then the orig.ident
  filename <- 'filtered_feature_bc_matrix.h5'
  if(!filter.spots) {
    filename <- 'raw_feature_bc_matrix.h5'
  }
  obj <- Seurat::Load10X_Spatial(spaceranger.dir, filename = filename, filter.matrix = filter.spots, slice = sample.name)
  obj$orig.ident <- sample.name
  tissue.positions <- get.tissue.position.metadata(spaceranger.dir)
  tissue.positions$spot_type <- "background"
  tissue.positions[tissue.positions$tissue == 1, "spot_type"] <- "tissue"
  obj <- AddMetaData(obj, tissue.positions)
  obj  
}

#' Create Seurat objects, one for eachVisium 10X sample / dataset.
#'  
#' @param spaceranger.dirs A named list, where each entry corresponds to a sample and provides 
#'                         a string directory holding the files filtered_feature_bc_matrix.h5 and
#'                         raw_feature_bc_matrix.h5 and the subdirectory spatial. This is likely
#'                         the outs directory created by spaceranger count or the equivalent.
#' @param filter.spots A boolean indicating whether to include all spots (filter.spots = FALSE) or only
#'                     those spots overlapping the tissue (filter.spots = TRUE)
#' @return A named list of Seurat objects, each annotated with the position and status 
#'         (spot_type = "tissue" or "background") for each spot. The names are those of spaceranger.dirs.
create.visium.seurat.objects <- function(spaceranger.dirs, filter.spots = TRUE) {
  samples <- names(spaceranger.dirs)
  names(samples) <- samples
  objs <-
    llply(samples, .fun = function(sample.name) {
      create.visium.seurat.object(sample.name, spaceranger.dirs[[sample.name]], filter.spots = filter.spots)
    })
  objs
}

#' Extract per-spot number/fraction of intronic, exonic, and intergenic (mapped) reads and number unmapped reads.
#' 
#' @param bam.file A string providing the path to an (indexed) bam file
#' @return A data.frame in which each row corresponds to a spot.
#'  It has rownames corresponding to the barcode of the spot and columns:
#'  CB: The sequence of the barcode associated to the spot.
#'  E: The number of exonic reads in the spot
#'  I: The number of intronic reads in the spot
#'  N: The number of intergenic reads in the spot
#'  tot.mapped.reads: The number of mapped reads in the spot
#'  tot.unmapped.reads: The number of unmapped reads in the spot
#'  tot.reads: tot.mapped.reads + tot.unmapped.reas
#'  E.frac.mapped: E / tot.mapped.reads
#'  I.frac.mapped: I / tot.mapped.reads
#'  N.frac.mapped: N / tot.mapped.reads
get.per.spot.alignment.metrics <- function(bam.file) {
  # The following tries to efficiently parse the bam by chromosome,
  # but this doesn't appear to work for FFPE/targeted sequencing.
  # For some reason, the rname/chromosome is NA for the vast majority (~98%)
  # of reads. Probably they are mapped to probe contigs/"chromosomes", which
  # aren't reflected in the rname.
  if(FALSE) {
    # Get the chromosome names
    ret <- scanBamHeader(c(bam.file))
    # The two elements returned by scanBamHeader are targets and text.
    # According to the Rsamtools docs:
    # The targets element contains target (reference) sequence lengths. The text element is
    # itself a list with each element a list corresponding to tags (e.g., ‘@SQ’) found in the header, and the
    # associated tag values.
    #sq.indices <- names(ret[[1]]$text) == "@SQ"
    #seqnames <- as.vector(unlist(lapply(ret[[1]]$text[sq.indices], function(x) gsub(x[1], pattern="SN:", replacement=""))))
    seq.lengths <- as.vector(ret[[1]]$targets)
    names(seq.lengths) <- names(ret[[1]]$targets)
    seq.names <- names(seq.lengths)
    names(seq.names) <- seq.names
    
    # For computational/memory efficiency, iterate over the chromosomes as opposed
    # to processing the entire file
    res <- ldply(seq.names, .parallel = FALSE,
                 .fun = function(seqname) {
                   print(seqname)
                   # CB: error-corrected spot barcode
                   # UB: error-corrected molecular barcode
                   # RE: E = exonic, N = intronic, I = intergenic
                   param <- ScanBamParam(tag=c('CB','RE'), which=GRanges(seqname, IRanges(1, seq.lengths[[seqname]])), flag=scanBamFlag(isUnmappedQuery=FALSE))
                   x <- scanBam(bam.file, param=param)[[1]]
                   if(is.null(x$tag$CB)) { return(NULL) }
                   if(is.null(x$tag$RE)) { return(NULL) }
                   df <- as.data.frame(x$tag)
                   # Combine counts within spot (for this chromosome)
                   ddply(df, .variables = c("CB"),
                         .fun = function(df.cb) {
                           if(nrow(df.cb) == 0) { return(NULL) }
                           data.frame(E = length(which(!is.na(df.cb$RE) & (df.cb$RE == "E"))),
                                      I = length(which(!is.na(df.cb$RE) & (df.cb$RE == "I"))),
                                      N = length(which(!is.na(df.cb$RE) & (df.cb$RE == "N"))),
                                      tot.mapped.reads = nrow(df.cb))
                         })
                 })
    
  } # end if(FALSE)
  param <- ScanBamParam(tag=c('CB','RE'), flag=scanBamFlag(isUnmappedQuery=FALSE))
  x <- scanBam(bam.file, param=param)[[1]]
  df <- as.data.frame(x$tag)
  # Combine counts within spot (for this chromosome)
  res <- ddply(df, .variables = c("CB"),
               .fun = function(df.cb) {
                 if(nrow(df.cb) == 0) { return(NULL) }
                 data.frame(E = length(which(!is.na(df.cb$RE) & (df.cb$RE == "E"))),
                            I = length(which(!is.na(df.cb$RE) & (df.cb$RE == "I"))),
                            N = length(which(!is.na(df.cb$RE) & (df.cb$RE == "N"))),
                            tot.mapped.reads = nrow(df.cb))
               })
  
  # Combine results for a spot (CB) across chromosomes
  spot.res <- ddply(res, .variables=c("CB"), .fun = function(df) data.frame(E = sum(df$E), I = sum(df$I), N = sum(df$N), tot.mapped.reads = sum(df$tot.mapped.reads)))
  
  # See description of _cell_ranger filtering here:
  # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview
  # I believe the bam we are reading already has duplicate UMIs filtered
  param <- ScanBamParam(tag=c('CB'),flag=scanBamFlag(isUnmappedQuery=TRUE))
  unmr <- scanBam(bam.file, param=param)[[1]]
  unmr.tbl <- as.data.frame(table(na.omit(unmr$tag$CB)))
  colnames(unmr.tbl) <- c("CB", "tot.unmapped.reads")
  
  spot.res <- merge(spot.res, unmr.tbl, all.x=TRUE)
  #spot.res <- spot.res[!is.na(spot.res$CB),]
  spot.res$tot.reads <- spot.res$tot.mapped.reads + spot.res$tot.unmapped.reads
  spot.res$E.frac.mapped <- spot.res$E / spot.res$tot.mapped.reads
  spot.res$I.frac.mapped <- spot.res$I / spot.res$tot.mapped.reads
  spot.res$N.frac.mapped <- spot.res$N / spot.res$tot.mapped.reads
  # tmp <- merge(unfiltered.objs[[3]][[]], spot.res, by.x = c("row.names"), by.y = c("CB"))
  # spot.res[!is.na(spot.res$CB),"CB"] <- "ambiguous.spot"
  # rownames(spot.res) <- spot.res$CB
  # spot.res[, !(colnames(spot.res) == "CB")]
  spot.res
}

#' Extract per-spot number of UMIs, total number of reads, and saturation.
#' 
#' Number of UMIs should match nCount_Spatial in metadata.
#' 
#' Saturation = PCR duplication = # non-unique reads / tot # reads
#' 
#' @param mol.info.file A string providing the path to "molecule_info.h5" file
#' @return A data.frame in which each row corresponds to a spot.
#'  It has rownames corresponding to the barcode of the spot and columns:
#'  num.reads: total number of reads to barcode
#'  num.umis: number of unique reads (UMIs) to barcode
#'  num.nonzero.genes: number of genes to barcode with one or more reads
#'  saturation: (num.reads - num.unis) / num.reads
get.per.spot.saturation.stats <- function(mol.info.file, filter = TRUE) {
  # Get the barcode and count of each UMI (ignoring the feature/gene)
  df <- data.frame(barcode = h5read(mol.info.file, "/barcodes")[h5read(mol.info.file, "/barcode_idx")+1], 
                   count = h5read(mol.info.file, "/count"),
                   feature.idx = h5read(mol.info.file, "/feature_idx"))
  # This is how you would extract the associated gene
  # gene <- h5read(mol.info.file, "/features/name")[h5read(mol.info.file, "/feature_idx")+1]
  # This should tell us whether the read mapped to an intron or exon, but I only see exonic reads
  # umi.type <- h5read(mol.info.file, "/umi_type")
  
  ls.info <- h5ls(mol.info.file)
  probe.set.key.flag <- grepl(ls.info$name, pattern="Probe Set")
  probe.set <- NULL
  if(length(which(probe.set.key.flag)) == 1) {
    probe.set <- h5read(mol.info.file, paste0(ls.info[probe.set.key.flag,"group"], "/", ls.info[probe.set.key.flag,"name"]))
  }
  res <- ddply(df, 
               .variables = c("barcode"), 
               .fun = function(tbl) {
                 # As stated here:
                 # https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
                 # non-targeted genes are removed from the 'filtered' (but not the 'unfiltered') matrix
                 if(filter && !is.null(probe.set)) {
                   tbl <- subset(tbl, feature.idx %in% probe.set)
                 }
                 num.reads <- sum(tbl$count)
                 num.umis <- nrow(tbl)
                 num.nonzero.genes <- length(unique(tbl$feature.idx[tbl$count > 0]))
                 saturation <- (num.reads - num.umis) / num.reads
                 data.frame(num.reads = num.reads, num.umis = num.umis, num.nonzero.genes = num.nonzero.genes, saturation = saturation)
               })
  rownames(res) <- res$barcode
  res <- res[, !(colnames(res) == "barcode")]
  res
}

#' Identify biotype of each gene in a a geneset
#' 
#' @param gene_set A character array carrying the querried gene symbols names
#' @param gene_db A Mart object database
#' @return A data.frame giving the biotype of each gene in the gene_set
get.biotypes_ <- function(gene_set,gene_db){
  gb <- getBM(attributes=c("external_gene_name", "gene_biotype", "chromosome_name"),filters = c("external_gene_name"), values=gene_set, mart=gene_db)
  # Label mitochondrial genes as any encoded on the "MT" chromosome
  flag <- grepl(gb$chromosome_name, pattern="MT", ignore.case = TRUE)
  gb[flag, "gene_biotype"] <- "MT"
  # Label ribosomal proteins (RP) as any with gene name Rps or Rpl.
  # Note that these are not the same as ribosomal (t)RNAs, which should be listed as rRNA by biomart
  flag <- grepl(gb$external_gene_name, pattern="^Rp[sl]", ignore.case = TRUE) # for mouse
  flag <- flag | grepl(gb$external_gene_name, pattern="^RP[SL]", ignore.case = TRUE) # for human
  gb[flag, "gene_biotype"] <- "RP"
  colnames(gb) <- c("gene", "gene_biotype", "chromosome_name")
  gb <- gb[, c("gene", "gene_biotype")]
  df <- data.frame(gene = gene_set)
  df <- merge(df, gb, all.x = TRUE)
  df
}

#' Simplify biotypes by condensing several into large categories
#' 
#' @param df A data.frame with a biotype column (e.g., with each row corresponding to a gene)
#' @param biotype.col The name of the biotype column in df
#' @return df modified so that immunoglobulin-related biotypes (IG_) are condensed into an IG biotype,
#'            T-cell receptor-related biotypes (TR_) are condensed into a TR biotype,
#'            and NA biotypes are converted to unknown.
condense.biotypes <- function(df, biotype.col = "gene_biotype") {
  flag <- is.na(df[, biotype.col])
  df[flag, biotype.col] <- "unknown"
  flag <- grepl(df[, biotype.col], pattern="IG_")
  df[flag, biotype.col] <- "IG"
  flag <- grepl(df[, biotype.col], pattern="TR_")
  df[flag, biotype.col] <- "TR"
  flag <- grepl(df[, biotype.col], pattern="pseudogene")
  df[flag, biotype.col] <- "pseudogene"
  df
}

#' Get the position information for each spot in a tissue
#' 
#' A function that extracts spot position information from the spaceranger output corresponding to 
#' a single tissue. See https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images
#' for a full description.
#'
#' @param spaceranger_dir The path to the spaceranger output.
#' @return A data.frame with columns:
#'  barcode: The sequence of the barcode associated to the spot.
#'  in_tissue: Binary, indicating if the spot falls inside (1) or outside (0) of tissue.
#'  array_row: The row coordinate of the spot in the array from 0 to 77. The array has 78 rows.
#'  array_col: The column coordinate of the spot in the array. In order to express the orange crate arrangement of the spots, this column index uses even numbers from 0 to 126 for even rows, and odd numbers from 1 to 127 for odd rows. Notice then that each row (even or odd) has 64 spots.
#'  pxl_row_in_fullres: The row pixel coordinate of the center of the spot in the full resolution image.
#'  pxl_col_in_fullres: The column pixel coordinate of the center of the spot in the full resolution image.
get.tissue.position.metadata <- function(spaceranger_dir) {
  image.dir <- paste0(spaceranger_dir, "/", "spatial/")
  tissue.positions <- 
    read.csv(file = file.path(image.dir, "tissue_positions_list.csv"), 
             col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, 
             as.is = TRUE, row.names = 1)
  tissue.positions
}