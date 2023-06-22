# Modified version of GetFit from
# https://raw.githubusercontent.com/saketkc/scRNA_NB_comparison/master/code/01_run_glm.R")
# to add parallelism (via llply)
my.GetFit <- function(cm, ncells = dim(cm)[2], type = "poisson") {
  cm <- cm[rowSums(cm) > 0, colSums(cm) > 0]
  ncells <- min(ncells, dim(cm)[2])
  total_umi <- colSums(cm)
  if (ncells == dim(cm)[2]) {
    # No sampling
    sampled_cells <- colnames(cm)
  } else {
    sampled_cells <- dds(total_umi, colnames(cm), ncells)
  }
  
  cm_sampled <- cm[, sampled_cells]
  total_umi_sampled <- colSums(cm_sampled)
  
  genes <- rownames(cm_sampled)
  names(genes) <- genes
  #pvals_list <- mclapply(genes, FUN = function(gene) {
  #  gene_umi <- as.vector(cm_sampled[gene, ])
  
  # pvals <- RunGLM(gene_umi = gene_umi, cell_umi = total_umi_sampled, type = type)
  # pvals
  #}, mc.cores = mc.cores)
  #  pvals_list <- lapply(genes, FUN = 
  pvals_list <- llply(genes, .parallel = TRUE, .fun =
                        function(gene) {
                          gene_umi <- as.vector(cm_sampled[gene, ])
                          pvals <- RunGLM(gene_umi = gene_umi, cell_umi = total_umi_sampled, type = type)
                          pvals
                        })
  names(pvals_list) <- genes
  
  pvalsfit <- bind_rows(pvals_list, .id = "gene")
  meanvar <- MeanVarFit(cm)
  meanvar$total_gene_umi <- rowSums(cm)
  
  meanvar$total_cell_umi <- sum(total_umi_sampled)
  meanvar$median_cell_umi <- median(total_umi_sampled)
  meanvar$mean_cell_umi <- mean(total_umi_sampled)
  
  meanvarfit <- inner_join(meanvar, pvalsfit, by = "gene")
  
  qval_cols <- c("qval_ssr_deviance", "qval_ssr_pearson", "qval_ssr_quantile", "qval_ssrstd_deviance", "qval_ssrstd_pearson", "qval_ssrstd_quantile")
  for (qval_col in qval_cols) {
    pval_col <- gsub("qval", "pval", qval_col)
    meanvarfit[, qval_col] <- GetQvalues(meanvarfit[, pval_col])
  }
  
  
  return(list(meanvarfit = meanvarfit, cell_attr = CellSummary(cm_sampled), gene_attr = GeneSummary(cm_sampled)))
}
