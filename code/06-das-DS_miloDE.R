suppressPackageStartupMessages({
  library(miloR)
  library(miloDE)
  library(SingleCellExperiment)
})

fun <- \(x) {
  y <- assign_neighbourhoods(x, reducedDim_name="PCA")
  res <- de_test_neighbourhoods(y, 
      verbose=FALSE, 
      design=~group_id, 
      sample_id="sample_id", 
      covariates="group_id") |> na.omit()
  if (!is.null(res)) {
    old <- c("pval", "pval_corrected_across_genes", "gene")
    new <- c("p_val", "p_adj", "gene_id")
    idx <- match(old, names(res))
    names(res)[idx] <- new
    is <- apply(z <- nhoods(y), 2, \(.) seq_len(ncol(x))[. == 1])
    fq <- sapply(is, \(.) max(prop.table(table(x[, .]$group_id))))
    res$n_cells <- colSums(z)[res$Nhood_center]
    res$n_nhood <- length(unique(res$Nhood))
    res$p_nhood <- fq[res$Nhood_center]
    res$i_nhood <- res$Nhood
  }
  return(res)
}