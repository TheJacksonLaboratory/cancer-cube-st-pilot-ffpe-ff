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
#' @return A data frame whose rows are spots, whose columns are deconvolved populations, and whose 
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    Also adds columns x and y, holding spatial coordinates of spots.
format.rctd.output_ <- function(rctd) {
  barcodes <- colnames(rctd@spatialRNA@counts)
  weights <- rctd@results$weights
  norm_weights <- normalize_weights(weights)
  df <- as.data.frame(norm_weights)
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
#' @return A data.frame containing the estimated quantiles for each probability value in quantiles (one row, and as many columns as values in quantiles)
calculate.expression.quantiles <- function(expr.mat, quantiles = seq(0, 1, by=0.1)) {
  all.means <- rowMeans(expr.mat)
  quantile(all.means, probs=quantiles)
}