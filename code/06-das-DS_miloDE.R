suppressPackageStartupMessages({
  library(miloDE)
  library(SingleCellExperiment)
  library(miloR)
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
    res$cell_id <- sapply(res$Nhood, \(i) list(which(nhoods(y)[,i]==1)))
    res$prop <- sapply(res$cell_id, 
      \(i) max(prop.table(table(colData(x)[i, "group_id"]))))
    res$nCells <- sapply(res$cell_id, length)
    res$nSubpop <- length(unique(res$Nhood))
    res$cell_id <- NULL
    res
  }
  return(res)
}