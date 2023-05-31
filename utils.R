
#' Run RCTD deconvolution by first ensuring that single-cell RNA-seq reference only includes genes in ST data.
#'
#' @param st.obj A Seurat object for the ST data
#' @param sc.counts The reference scRNA-seq count matrix (with genes as rows and cells as columns). 
#'                  Genes should be in the same namespace (e.g., symbols or Ensembl ids) as the ST count data in st.obj.
#' @param sc.cell.types The cell type assignments of each cell in sc.counts stored in a vector, named according to the columns/cells of sc.counts
#' @param rds.output.file If non-null, an output file name to store the RCTD results. Default: NULL.
#'                        If the file exists, the output will be read and returned by this function.
#' @param run.rctd Boolean indicating whether to actually deconvolve data (run.RCTD) or simply to create.RCTD.
#' @param intersect.sc.and.st.genes Boolean indicating whether to subset scRNA-seq genes to those within
#'        ST data.
#' @return An RCTD object
#' 
#' Assuming st.obj is the filtered feature obj, the associated count matrix should _already_ be filtered to include only
#' targeted genes, as described here:
#' https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
rctd.wrapper <- function(st.obj, sc.counts, sc.cell.types, rds.output.file = NULL, run.rctd = TRUE, intersect.sc.and.st.genes = FALSE) {
  if(!is.null(rds.output.file) && file.exists(rds.output.file)) {
    myRCTD <- readRDS(rds.output.file)
    gc()
    return(myRCTD)
  }
  
  # Get raw ST coiunts
  st.counts <- GetAssayData(st.obj, assay="Spatial", slot="counts")
  
  # Get the spot coordinates
  st.coords <- st.obj[[]][, c("col", "row")]
  colnames(st.coords) <- c("x","y")
  
  # Get the # of UMI per spot
  nUMI <- colSums(st.counts) # In this case, total counts per pixel is nUMI
  
  # Create the RCTD 'puck', representing the ST data
  puck <- SpatialRNA(st.coords, st.counts, nUMI)
              
  # Create an RCTD scRNA-seq reference. 
  # NB: we do this for each of potentially multiple ST samples because some of those
  # samples may be from FFPE assays and others from fresh frozen. These will have
  # different gene sets that should be reflected (i.e., filtered) in the scRNA-seq reference (sc.counts).
  
  if(intersect.sc.and.st.genes) {
    stop("Would not recommend restricting scRNA-seq genes to those in ST!")
    # spacexr/RCTD:get_de_genes has an expression threshold (expr_threshold)
    # that any marker gene must exceed. Restricting the set of genes (to those in ST),
    # impacts their normalization (by total UMI) in RCTD. I have observed that this
    # can lead to significant differences in the number of genes passing the 
    # expression threshold (e.g., twice as many following subsetting the genes).
    # Note that get_de_genes will only consider as markers those genes in common
    # anyways between scRNA-seq and ST. As such, it intuitively makes sense
    # to consistently define the normalized expression values of the scRNA-seq data
    # irrespective of the ST data -- i.e., without subsetting.
    sc.counts <- sc.counts[rownames(sc.counts) %in% rownames(st.counts),]
  }
  sc.cell.types <- sc.cell.types[colnames(sc.counts)]
  nUMI <- colSums(sc.counts)
  reference <- Reference(sc.counts, sc.cell.types, nUMI, n_max_cells = ncol(sc.counts) + 1)
              
  max_cores <- min(5, detectCores() - 1)
  
  # Note that we have keep_reference = TRUE here. Without it, the default is
  # keep_reference = FALSE, which would ignore the Reference we have 
  # intentionally created above.
  myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, keep_reference = TRUE)
  
  # NB: running in full mode!
  if(run.rctd) {
    myRCTD <- suppressPackageStartupMessages(run.RCTD(myRCTD, doublet_mode = 'full'))
  }
  
  if(!is.null(rds.output.file)) {
    saveRDS(myRCTD, rds.output.file)
  }
  
  gc()
  myRCTD
}

#' Create a .csv file for the automatic analysis of ST data, as it is described 
#' here: https://github.com/sdomanskyi/spatialtranscriptomics
#' 
#' @param base_dir A base directory where the .csv file will be saved
#' @param species_sample Vector describing the species of each sample
#' 
#' 
#' # Write csv file for automatic ST pipeline 
write.automated.csv.file <- function(base_dir, dataset_names, species_sample) {
  temp_var<- paste0(base_dir,dataset_names,'/spaceranger/')
  temp_var2<- replicate(length(dataset_names), "")
  ST_automated<- data.frame (sample_id= names(dataset_names), species= species_sample, st_data_dir=temp_var, sc_data_dir= temp_var2 )
  dir.create(file.path(base_dir, '/automatic_analysis/'), recursive=TRUE)
  temp <- temp <- strsplit(base_dir,"/")
  curr_datatypes <- tail(temp[[1]], n=1)
  write.csv(ST_automated,paste0(base_dir,'/automatic_analysis/',curr_datatypes,'.csv'), row.names = FALSE ,quote=FALSE)
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
  
  file <- file.path(image.dir, "tissue_positions_list.csv")
  col.names <- c("barcodes", "tissue", "row", "col", "imagerow", "imagecol")
  if(file.exists(file)) {
    tissue.positions <- 
      read.csv(file = file, col.names = col.names, header = FALSE, as.is = TRUE, row.names = 1)
    return(tissue.positions)
  }
  file <- file.path(image.dir, "tissue_positions.csv")
  if(!file.exists(file)) {
    stop("Found neither tissue_positions_list.csv nor tissue_positions.csv")
  }
  tissue.positions <- 
    read.csv(file = file, header = TRUE, as.is = TRUE)
  colnames(tissue.positions) <- col.names
  return(tissue.positions)
}

#' Create a Seurat object based on the spaceranger baseline file
#'
#' @param spaceranger_dirs The directory where the folders for each dataset exist
#' @param dataset The name of queried dataset
#' @param spaceranger_bsl_filename The spaceranger file where the details for each Seurat data odject are extracted
#' @param filter_param Binary variable for filtering tissue positions
#' @return A Suerat object
create.Seurat.object <- function(spaceranger_dirs, dataset, spaceranger_bsl_filename, filter_param){
  obj <- Seurat::Load10X_Spatial(spaceranger_dirs[[dataset]], filename = spaceranger_bsl_filename, filter.matrix = filter_param, slice = dataset)
  obj$orig.ident <- dataset
  tissue.positions <- get.tissue.position.metadata(spaceranger_dirs[[dataset]])
  tissue.positions$spot_type <- "background"
  tissue.positions[tissue.positions$tissue == 1, "spot_type"] <- "tissue"
  obj <- AddMetaData(obj, tissue.positions)
  obj
}

#' Merge all Seurat objects in one object using harmony (https://portals.broadinstitute.org/harmony/articles/quickstart.html)
#' and cluster the data based on gene expression data and RUN UMAP
#' 
#' @param objs A merged Seurat object including all Seurat objects we want to merge
#' @param norm_param Choose which normalization assay will be used
#' @return merged_obj One Seurat object containing all unique Seurat objs
merge.Seurat.objs <- function(objs, norm_param){
  merged_obj <- ScaleData(merged_obj, assay=norm_param)
  VariableFeatures(merged_obj, assay=norm_param) <- unique(unname(unlist(lapply(filtered.objs, function(x) VariableFeatures(x, assay=norm_param)))))
  merged_obj <- RunPCA(merged_obj, verbose = FALSE)
  ElbowPlot(merged_obj, ndims = 40)
  mat <- Seurat::GetAssayData(merged_obj, assay = norm_param)
  pca <- merged_obj[["pca"]]
  
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(as.matrix(mat)))
  
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  cumVarExplained <- cumsum(varExplained)
  
  merged_obj <- RunHarmony(merged_obj, group.by.vars = make.names("orig ident"))
  
  # Pick the number of dimensions based on a flattening of the elbow plot
  merged_obj <- FindNeighbors(merged_obj,reduction = "harmony", dims = 1:20)
  merged_obj <- FindClusters(merged_obj, reduction = "harmony",verbose = FALSE)
  merged_obj <- RunUMAP(merged_obj, reduction = "harmony",dims = 1:20)
  merged_obj
}

#' Return count matrix from a Seurat object
#' 
#' @param obj A Seurat object.
#' @return A (sparse) matrix holding the spot counts (i.e., from the "counts" slot)
get.count.matrix <- function(obj) {
  as.matrix(x = GetAssayData(object = obj, assay = "Spatial", slot = "counts"))
}


#' Downsample each column of a matrix to a desired total count.
#'
#' A function that, independently for each column/cell, uses a multinomial distribution, parameterized 
#' by the frequency of each gene within that column/cell and the desired total count, to sample a vector of gene counts.
#' This is in contrast to related functions such as SampleUMI (in Seurat), which applies a binomial to each gene (!?)
#' and thus only approximates the desired total count, and downsampleMatrix (in DropletUtils), which downsamples to
#' a fraction of the original number of reads.
#'
#' @param expr_mat A (possibly sparse) matrix of counts.
#' @param tot_cnt The desired number of counts _for each column_
#' @return A downsample matrix
downsample.matrix <- function (expr_mat, tot_cnt) 
{
  down_sample <- function(x) {
    prob <- x/sum(x)
    return(rmultinom(1, tot_cnt, prob))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

#' Compute the number of expressed genes/features in a downsampled matrix.
#' 
#' A function to downsample the input matrix such that each cell/column has tot_cnt reads and
#' that reports the _median_ (across cells/columns) number of genes/features with non-zero counts.
#'
#' @param expr_mat A (possibly sparse) matrix of counts.
#' @param tot_cnt The desired number of counts _for each column_
#' @return The median (across cells/columns) number of genes/features with non-zero counts in the downsampled matrix.
compute.num.features.at.total.reads <- function(expr_mat, tot_cnt) {
  ds <- downsample.matrix(expr_mat, tot_cnt)
  return(median(colSums(ds > 0)))
}

#' Return most highly-expressed genes
#' 
#' A function to return the most highly-expressed genes (i.e., those with highest value) from a matrix.
#'
#' @param mat A (sparse) matrix.
#' @param n.top The number of genes to return.
#' @param assay The assay from which to access the gene exprssion/values in the Seurat object.
#' @param slot The slot in the assay form which to access the gene expression/values.
#' @return A vector of the highest-expressed genes.
get.top.genes.matrix <- function(mat, n.top = 20, assay = "Spatial", slot = "data") {
  rs <- rowSums(mat)
  rs <- rs[order(rs, decreasing=TRUE)]
  top.genes <- names(rs)[1:n.top]
  top.genes
}

#' Return most highly-expressed genes
#' 
#' A function to return the most highly-expressed genes (i.e., those with highest value) from a Seurat assay.
#'
#' @param obj A Seurat object.
#' @param n.top The number of genes to return.
#' @param assay The assay from which to access the gene exprssion/values in the Seurat object.
#' @param slot The slot in the assay form which to access the gene expression/values.
#' @return A vector of the highest-expressed genes.
get.top.genes <- function(obj, n.top = 20, assay = "Spatial", slot = "data") {
  mat <- GetAssayData(obj, assay = assay, slot = slot)
  get.top.genes.matrix(mat, n.top = n.top, assay = assay, slot = slot)
}

#' Return data subsetted by spot tissue status
#' 
#' A function to return the subset of the data corresponding to spots that are (or are not) in the tissue
#'
#' @param obj A Seurat object.
#' @param tissue.val Spots with tissue.val for the tissue.col metadata field will be retained.
#' @param tissue.col The metadata column of the Seurat object holding the tissue status
#' @param assay The assay from which to access the gene exprssion/values in the Seurat object.
#' @param slot The slot in the assay form which to access the gene expression/values.
#' @return A (sparse) matrix of spots that are (or are not) in tissue
subset.matrix.based.on.tissue.status <- function(obj, tissue.val, tissue.col = "tissue", assay = "Spatial", slot = "data") {
  cell.subset <- rownames((obj@meta.data)[(obj@meta.data)[,tissue.col] == tissue.val,])
  mat <- GetAssayData(obj, assay = assay, slot = slot)
  mat[,cell.subset]
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

#' Identify frequency of each biotype within a geneset
#' 
#' @param gene_set A character array carrying the querried gene symbols names
#' @param gene_db A Mart object database
#' @return A data.frame giving the count (Freq) of each biotype in gene_set
get.biotypes <- function(gene_set,gene_db){
  gb <- get.biotypes_(gene_set, gene_db)
  querried.biotypes <- as.data.frame(table(gb$gene_biotype))
  colnames(querried.biotypes) <- c("biotype", "Freq")
  o <- order(querried.biotypes$Freq, decreasing=TRUE)
  querried.biotypes <- querried.biotypes[o,]
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

#' Try to load package, installing first (via pacman) if necessary
#' 
#' @param package A string representing the package to load
#' @return Nothing
p_load_ <- function(package) {
  if(!require(pacman)) {
    install.packages(pacman)
  }
  suppressPackageStartupMessages(p_load(package))
}

#' Create a panel of plots, each showing deconvolved fractions of a particular population.
#' 
#' @param rctd An RCTD object, from the spacexr package.
#' @param normalize Boolean indicating whether weights should be normalized to sum to 1.
#' @return A data frame whose rows are spots, whose columns are deconvolved populations, and whose 
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    Also adds columns x and y, holding spatial coordinates of spots.
format.rctd.output_ <- function(rctd, normalize = TRUE) {
  barcodes <- colnames(rctd@spatialRNA@counts)
  weights <- rctd@results$weights
  if(normalize) {
    weights <- normalize_weights(weights)
  }
  df <- as.data.frame(weights)
  df$x <- rctd@spatialRNA@coords$x
  df$y <- rctd@spatialRNA@coords$y
  df
}

#' Parse ligand-receptor interactions from cellphonedb.
#' 
#' @return A list with the following entries:
#'         interaction.tbl: a data.frame where rows are ligand-receptor interactions and with the following columns:
#'                          id_cp_interaction: cellphonedb name for ligand-receptor interaction pair
#'                          partner_a: cellphonedb name for partner a in the interaction (i.e., ligand)
#'                          partner_b: cellphonedb name for partner b in the interaction (i.e., receptor)
#'                          ligand_hgnc_symbols: comma-separated list of ligand hgnc symbols
#'                          ligand_ensembl_ids: comma-separated list of ligand ensembl ids
#'                          receptor_hgnc_symbols: comma-separated list of ligand hgnc symbols
#'                          receptor_ensembl_ids: comma-separated list of ligand ensembl ids
#'         all.ligand.hgnc.symbols: a vector of all ligand hgnc symbols (i.e., concatenation of ligand_hgnc_symbols from interaction.tbl)  
#'         all.ligand.ensembl.ids: a vector of all ligand ensembl ids (i.e., concatenation of ligand_ensembl_ids from interaction.tbl)  
#'         all.receptor.hgnc.symbols: a vector of all receptor hgnc symbols (i.e., concatenation of receptor_hgnc_symbols from interaction.tbl)  
#'         all.receptor.ensembl.ids: a vector of all receptor ensembl ids (i.e., concatenation of receptor_ensembl_ids from interaction.tbl)  
read.cellphonedb.ligand.receptors <- function() {
  # Read the cellphonedb tables
  interaction.tbl <- read.table("https://raw.githubusercontent.com/ventolab/CellphoneDB/master/cellphonedb/src/core/data/interaction_input.csv", sep=",", header=TRUE)
  complex.tbl <- read.table("https://raw.githubusercontent.com/ventolab/CellphoneDB/master/cellphonedb/src/core/data/complex_input.csv", sep=",", header=TRUE)
  gene.tbl <- read.table("https://raw.githubusercontent.com/ventolab/CellphoneDB/master/cellphonedb/src/core/data/gene_input.csv", sep=",", header=TRUE)
  protein.tbl <- read.table("https://raw.githubusercontent.com/ventolab/CellphoneDB/master/cellphonedb/src/core/data/protein_input.csv", sep=",", header=TRUE)
  
  # The complexes table has a complex_name (e.g., ACVL1_BMPR2) and up to four associated proteins 
  # (identified by uniprot ids uniprot_{1,4}).
  # To translate these proteins to genes, merge with the genes table, which has columns
  # gene_name, uniport, hgnc_symbol, and ensembl
  uniprot.cols <- grep(x = colnames(complex.tbl), pattern="uniprot", value = TRUE)
  suffixes <- gsub(uniprot.cols, pattern = "uniprot", replacement="")
  
  for(i in 1:length(uniprot.cols)) {
    gene.i.tbl <- gene.tbl
    colnames(gene.i.tbl) <- paste0(colnames(gene.i.tbl), suffixes[i])
    complex.tbl <- merge(complex.tbl, gene.i.tbl, all.x = TRUE)
  }
  
  # Concatenate all of the hgnc symbol and ensembl columns
  hgnc.symbol.cols <- grep(x = colnames(complex.tbl), pattern="hgnc_symbol", value = TRUE)
  complex.tbl[, "hgnc_symbols"] <- unlist(apply(complex.tbl[, hgnc.symbol.cols], 1, function(row) paste0(na.omit(row), collapse=",")))
  ensembl.cols <- grep(x = colnames(complex.tbl), pattern="ensembl", value = TRUE)
  complex.tbl[, "ensembl_ids"] <- unlist(apply(complex.tbl[, ensembl.cols], 1, function(row) paste0(na.omit(row), collapse=",")))
  
  # Now merge the interaction and complex tables -- where partner_{a,b} in the interaction table is the {ligand,receptor} complex (or uniprot protein id)
  suffixes <- c("_a", "_b")
  for(i in 1:length(suffixes)) {
    complex.i.tbl <- complex.tbl[, c("complex_name", "hgnc_symbols", "ensembl_ids")]
    colnames(complex.i.tbl) <- paste0(colnames(complex.i.tbl), "_complex", suffixes[i])
    interaction.tbl <- merge(interaction.tbl, complex.i.tbl, all.x = TRUE, by.x = paste0("partner", suffixes[i]), by.y = paste0("complex_name_complex", suffixes[i]))
  }
  
  # Now merge the interaction and protein tables -- to capture cases in which partner_{a,b} is a uniprot id
  # But, first merge the protein and gene tables
  protein.tbl <- merge(protein.tbl, gene.tbl, all.x = TRUE, by = "uniprot")
  for(i in 1:length(suffixes)) {
    protein.i.tbl <- protein.tbl[, c("uniprot", "hgnc_symbol", "ensembl")]
    colnames(protein.i.tbl) <- paste0(colnames(protein.i.tbl), "_protein", suffixes[i])
    interaction.tbl <- merge(interaction.tbl, protein.i.tbl, all.x = TRUE, by.x = paste0("partner", suffixes[i]), by.y = paste0("uniprot_protein", suffixes[i]))
  }
  
  # At this point, we have potentially multiple rows per interaction and multiple columns with both hgnc symbols and ensembl ids.
  # Let's concatentate all of these
  concat.vals.across.rows.and.cols <- function(df, pattern) {
    cols <- grep(colnames(df), pattern=pattern, value=TRUE)
    vals <- unlist(apply(df[, cols, drop=FALSE], 1, function(row) paste0(na.omit(row), collapse=",")))
    vals <- paste0(vals, collapse=",")
    vals <- sort(unique(unlist(strsplit(vals, split=",")[[1]])))
    vals <- paste0(vals, collapse=",")
    vals
  }
  
  interaction.tbl <-
    ddply(interaction.tbl, .variables = c("id_cp_interaction", "partner_a", "partner_b"),
          .fun = function(df) {
            lig.syms <- concat.vals.across.rows.and.cols(df, pattern = ".*symbol.*_a")
            rec.syms <- concat.vals.across.rows.and.cols(df, pattern = ".*symbol.*_b")
            lig.ids <- concat.vals.across.rows.and.cols(df, pattern = ".*ensembl.*_a")
            rec.ids <- concat.vals.across.rows.and.cols(df, pattern = ".*ensembl.*_b")
            data.frame("ligand_hgnc_symbols" = lig.syms, "ligand_ensembl_ids" = lig.ids,
                       "receptor_hgnc_symbols" = rec.syms, "receptor_ensembl_ids" = rec.ids)
          })
  
  all.ligand.symbols <- strsplit(concat.vals.across.rows.and.cols(interaction.tbl, pattern="ligand_hgnc_symbols"), split=",")[[1]]
  all.ligand.ids <- strsplit(concat.vals.across.rows.and.cols(interaction.tbl, pattern="ligand_ensembl_ids"), split=",")[[1]]
  all.receptor.symbols <- strsplit(concat.vals.across.rows.and.cols(interaction.tbl, pattern="receptor_hgnc_symbols"), split=",")[[1]]
  all.receptor.ids <- strsplit(concat.vals.across.rows.and.cols(interaction.tbl, pattern="receptor_ensembl_ids"), split=",")[[1]]
  
  lst <- list("interaction.tbl" = interaction.tbl,
              "all.ligand.hgnc.symbols" = all.ligand.symbols,
              "all.ligand.ensembl.ids" = all.ligand.ids,
              "all.receptor.hgnc.symbols" = all.receptor.symbols,
              "all.receptor.ensembl.ids" = all.receptor.ids
  )
  return(lst)
}

#' Calculate an "expression" for a ligand-receptor pair.
#' 
#' Represents the LR pair expression as the expression of a gene whose mean expression (across samples)
#' is minimal. This seems similar to the approach taken by cellphonedb, wherein the gene with minimal
#' expression within a multi-subunit heteromeric complex is used as the expression of that (ligand or receptor)
#' complex. Here, we are effectively taking the minimum of the ligand and receptor complexes for use as the
#' LR pair expression.
#' 
#' @param expr.mat An expression matrix (with no assumed units), whose rows are genes (in no particular namespace) and whose columns are samples/cells/spots.
#' @param ligand.genes A vector of ligand genes (in the same namespace as the rows of expr.mat)
#' @param receptor.genes A vector of receptor genes (in the same namespace as the rows of expr.mat)
#' @return A vector of expression for the ligand-receptor pair, named with the columns / samples of the expression matrix.
calculate.ligand.receptor.pair.expression <- function(expr.mat, ligand.genes, receptor.genes) {
  all.lr.genes <- unique(c(ligand.genes, receptor.genes))
  all.lr.genes <- all.lr.genes[all.lr.genes %in% rownames(expr.mat)]
  lr.means <- rowMeans(expr.mat[all.lr.genes,])
  # Find the gene with the minimum (average) expression. We will use this to represent the
  # expression of the lr pair
  gene.rep <- names(which.min(lr.means))[1]
  expr.mat[gene.rep,]
}

#' Calculate the quantiles of the mean expression (over samples/cell/spots/columns) of an expression matrix.
#'  
#' @param expr.mat An expression matrix (with no assumed units), whose rows are genes (in no particular namespace) and whose columns are samples/cells/spots.
#' @param quantiles A vector of probability values
#' @summary.func Function to apply to rows of expression matrix to summarize each gene's expression
#' @return A data.frame containing the estimated quantiles for each probability value in quantiles (one row, and as many columns as values in quantiles)
calculate.expression.quantiles <- function(expr.mat, quantiles = seq(0, 1, by=0.1), summary.func = mean) {
  # gene.summaries <- rowMeans(expr.mat)
  gene.summaries <- apply(expr.mat, 1, summary.func)
  quantile(gene.summaries, probs=quantiles)
}

#' Add metadata to Seurat object.
#'  
#' This function is a simple wrapper around AddMetaData
#' 
#' @param obj A Seurat object
#' @param metadata.df A data.frame expected to have all of the rownames (corresponding to cells/spots) in existing metadata for obj (obj[[]])
#' @return obj modified to hold the metadata in metadata.df
add.metadata.to.seurat.obj <- function(obj, metadata.df) {
  # rows <- rownames(obj[[]])
  # print(head(rows))
  # print(table(rows %in% rownames(metadata.df)))
  # stopifnot(rows %in% rownames(metadata.df))
  # metadata.df <- metadata.df[rows,]
  colnames(metadata.df) <- make.names(colnames(metadata.df))
  metadata.df <- metadata.df[, !(colnames(metadata.df) %in% colnames(obj[[]])), drop=FALSE]
  cols.to.add <- colnames(metadata.df)
  all.meta <- merge(obj[[]], metadata.df, by = "row.names", all.x = TRUE)
  rownames(all.meta) <- all.meta$Row.names
  ids <- Cells(obj)
  all.meta <- all.meta[ids,]
  obj <- AddMetaData(obj, all.meta[, cols.to.add, drop=FALSE])
  obj
}

#' Enable parallel execution.
#'  
#'  #' @param obj A Seurat object
#' @return The number of cores available for parallel execution.
setup.parallel.environment <- function(max.cores = NULL) {
  num.cores <- detectCores()
  if(!is.na(num.cores) && (num.cores > 1)) {
    if(!is.null(max.cores)) { num.cores <- min(num.cores, max.cores) }
    suppressPackageStartupMessages(p_load("doMC"))
    cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
    registerDoMC(cores=(num.cores-1))
  }
  num.cores
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
  # by Load10X_Spatial)
  
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

#' Concatenate Seurat objects into one object
#' 
#' @param objs A list of Seurat objects
#' @return A single Seurat object that is the concatenation of the input objects
concatenate.seurat.objects <- function(objs) {
  all.objs <- Reduce(merge, objs)
  Idents(all.objs) <- all.objs$orig.ident
  all.objs
}

#' Identify biotypes from a geneset given the orgnaism
#' 
#' @param gene_set A character array carrying the querried gene symbols names
#' @param gene_db A Mart object database
#' @return querried.biotypes
get.annotationhub.biotypes <- function(gene_set){
  stop("Not implemented")
  ah <- AnnotationHub()
  human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
  foo <- human_ens[["AH75011"]]
  df <- genes(foo, return.type = "data.frame")
  gb <- getBM(attributes=c("external_gene_name", "gene_biotype"),filters = c("external_gene_name"), values=gene_set, mart=gene_db)
  querried.biotypes <- as.data.frame(table(gb$gene_biotype))
  colnames(querried.biotypes) <- c("biotype", "Freq")
  o <- order(querried.biotypes$Freq, decreasing=TRUE)
  querried.biotypes <- querried.biotypes[o,]
}

#' Extract per-spot number/fraction of intronic, exonic, and intergenic (mapped) reads and number unmapped reads.
#' 
#' @param bam.file A string providing the path to an (indexed) bam file
#' @return A data.frame in which each row corresponds to a spot.
#'  It has rownames corresponding to the barcode of the spot and columns:
#'  CThe sequence of the barcode associated to the spot.
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

#' Extract per-spot number/fraction of intronic, exonic, and intergenic (mapped) reads and number unmapped reads over Seurat objects.
#' 
#' Wrapper to call get.per.spot.alignment.metrics across Seurat objects. 
#' 
#' @param spaceranger_dirs _Named_ list of spacernanger output directories, one per dataset.
#' @return A list of tables of alignment metrics, each output by get.per.spot.alignment.metrics
get.all.per.spot.alignment.metrics <- function(spaceranger_dirs, prefix = NULL) {
  nms <- names(spaceranger_dirs)
  names(nms) <- nms
  align.metrics <-
    llply(nms, .parallel = FALSE,
          .fun = function(nm) {
            print(nm)
            alignment.metric.file <- NULL
            if(!is.null(prefix)) {
              alignment.metric.file <- paste0(prefix, "/", nm, "-alignment-metrics.csv")
            }
            bam.file <- paste0(spaceranger_dirs[[nm]], "/possorted_genome_bam.bam")
            df <- NULL
            if(is.null(alignment.metric.file) || !file.exists(alignment.metric.file)) { 
              df <- get.per.spot.alignment.metrics(bam.file)
              write.table(file=alignment.metric.file, df, row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)
            } else {
              df <- as.data.frame(fread(alignment.metric.file))
            }
            df
          })
  align.metrics
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

#' Find the empirical distance between spots in the x and y directions.
#' 
#' It makes most sense for the input Seurat object to be unfiltered (by background vs tissue/foreground)
#' so that there aren't gaps between rows / columns. Nevertheless, the code attempts to
#' overcome this possibility by taking the most common distance between neighboring spots.
#' 
#' @param obj A Seurat object. 
#' @return A 2-tuple with elements:
#'           dx: the separation between neighboring spots in the x dimension
#'           dy: the separation between neighboring spots in the y dimension
get.visium.spot.distance.separation <- function(obj) {
  # Within each row, calculate the distance between spots ordered by column.
  tmp <- ddply(obj[[]], .variables=c("row"), 
               .fun = function(df) {
                 o <- order(df$col, decreasing=FALSE)
                 df <- df[o,]
                 diffs = unlist(lapply(1:(nrow(df)-1), function(i) df[i+1,"imagecol"] - df[i,"imagecol"]))
                 na.omit(data.frame(diff=diffs))
               })
  # Take the most common such distance between columns of neighboring spots as the
  # separation in the x direction -- dx
  tbl <- table(tmp$diff)
  dx <- as.numeric(names(sort(tbl, decreasing=TRUE))[1])
  # Within each col, calculate the distance between spots ordered by row.
  tmp <- ddply(obj[[]], .variables=c("col"), 
               .fun = function(df) {
                 o <- order(df$row, decreasing=FALSE)
                 df <- df[o,]
                 diffs = unlist(lapply(1:(nrow(df)-1), function(i) df[i+1,"imagerow"] - df[i,"imagerow"]))
                 na.omit(data.frame(diff=diffs))
               })
  # Take the most common such distance between rows of neighboring spots as the
  # separation in the y direction -- dy
  tbl <- table(tmp$diff)
  dy <- as.numeric(names(sort(tbl, decreasing=TRUE))[1])
  return(c(dx, dy))
}

#' Return the genes targeted by sequencing.
#' 
#' For fresh frozen, this will be all genes (in the Seurat object's expression matrix).
#' For FFPE, this will be those in the probe set.
#' 
#' @param mol.info.file A string providing the path to "molecule_info.h5" file
#' @return A list with named entries:
#'           gene.targets: A vector of gene names
#'           gene.off.targets: A vector of gene names
#'           gene.id.targets: A vector of gene ids
#'           gene.id.off.targets: A vector of gene ids
get.targeted.genes <- function(mol.info.file) {
  gene.ensg.ids <- h5read(mol.info.file, "/features/id")
  gene.names <- h5read(mol.info.file, '/features/name')
  
  ls.info <- h5ls(mol.info.file)
  probe.set.key.flag <- grepl(ls.info$name, pattern="Probe Set")
  if(length(which(probe.set.key.flag)) == 0) {
    # This is not a targeted assay
    return(list("gene.targets" = gene.names,
                "gene.off.targets" = NULL,
                "gene.id.targets" = gene.ensg.ids,
                "gene.id.off.targets" = NULL))
  }
  # Empirically, it appears that these indices are zero-based. I can't find that directly stated online, however:
  # 1. Other indices are zero-based, e.g., barcode_idx and feature_idx here:
  # https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/molecule_info
  # 2. As stated here
  # https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
  # the filtered feature for FFPE/targeted assays should only include the targeted genes. When I assume zero-based indices (as below),
  # I get a near perfect correspondence between probe.gene.names and the rownames of the (filtered) count matrix.
  # When I assume a one-based index, as I had previously, the correspondence is poor.
  probe.gene.indices <- h5read(mol.info.file, paste0(ls.info[probe.set.key.flag,"group"], "/", ls.info[probe.set.key.flag,"name"])) + 1
  gene.indices <- 1:length(gene.ensg.ids)
  nonprobe.gene.indices <- gene.indices[!(gene.indices %in% probe.gene.indices)]
  # probe.gene.indices <- gene.indices[(gene.indices %in% probe.gene.ids)]
  nonprobe.gene.names <- sort(gene.names[nonprobe.gene.indices])
  probe.gene.names <- gene.names[probe.gene.indices]
  nonprobe.gene.ids <- gene.ensg.ids[nonprobe.gene.indices]
  probe.gene.ids <- gene.ensg.ids[probe.gene.indices]
  list("gene.targets" = probe.gene.names,
       "gene.off.targets" = nonprobe.gene.names,
       "gene.id.targets" = probe.gene.ids,
       "gene.id.off.targets" = nonprobe.gene.ids)
}

#' Wrapper to apply sctransform to a seurat object.
#' 
#' @param obj A seurat object
#' @param use.v2 Boolean indicating whether to apply the "v2" flavor of sctransform
#' @param assay The assay to apply sctransform to
#' @return obj modified to hold sctransform'ed results in the "SCT" assay
apply.sctransform_ <- function(obj, use.v2 = TRUE, assay = "Spatial") {
  # obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  if(use.v2) { 
    obj <- SCTransform(obj, assay = assay, vars.to.regress = c(paste0("nCount_", assay)), verbose = FALSE, vst.flavor = "v2", method = "glmGamPoi")
  } else {
    obj <- SCTransform(obj, assay = assay, vars.to.regress = c(paste0("nCount_", assay)), verbose = FALSE)
  }
  return(obj)
}

#' Wrapper to apply sctransform to a list of seurat objects.
#' 
#' @param objs A list of seurat objects
#' @return the list of objs, with each entry modified to hold sctransform'ed results in the "SCT" assay
apply.sctransform <- function(objs, use.v2 = TRUE) {
  objs <-
    llply(objs,
          .fun = function(obj) {
            apply.sctransform_(obj, use.v2 = use.v2)
          })
  objs
}

