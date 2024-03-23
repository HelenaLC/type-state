suppressPackageStartupMessages({
  library(MouseGastrulationData)
  library(SingleCellExperiment)
  library(scater)
})

fun <- \() {
  # load SCE
  x <- WTChimeraData(samples=5:10)
  
  # standardize metadata
  colData(x) <- DataFrame(
    cluster_id=x$celltype.mapped,
    sample_id=x$sample,
    group_id=x$tomato)
  
  # drop doublets & undetected genes
  x <- x[, x$cluster_id != "Doublet"]
  x <- x[rowSums(counts(x)) > 0, ]
  
  # restrict to clusters with
  # >50 cells in >2 samples per group
  cs <- split(seq(ncol(x)), x$group_id)
  ks <- sapply(cs, \(.) {
    ns <- table(x$cluster_id[.], x$sample_id[.])
    names(which(rowSums(ns > 50) > 2))
  })
  ks <- Reduce(intersect, ks)
  dim(x <- x[, x$cluster_id %in% ks])
  table(x$cluster_id, x$sample_id)
  
  for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])
  
  # store gene/cell identifiers
  rowData(x)$gene_id <- rownames(x)
  colData(x)$cell_id <- colnames(x)
  
  return(x)
}
