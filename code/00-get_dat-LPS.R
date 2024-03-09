suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
})

fun <- \() {
  # load SCE
  x <- muscData::Crowell19_4vs4()
  
  # remove undetected genes
  x <- x[rowSums(counts(x)) > 0, ]
  
  # drop multiplet cells
  x <- scDblFinder(x, clusters="cluster_id", samples = "sample_id")
  x <- x[,x$scDblFinder.class == "singlet"]
  
  # break sample pairing
  colData(x) <- DataFrame(
    cluster_id=x$cluster_id,
    sample_id=x$sample_id,
    group_id=x$group_id)
  
  # restrict to clusters with
  # >30 cells in >2 samples per group
  cs <- split(seq(ncol(x)), x$group_id)
  ks <- sapply(cs, \(.) {
    ns <- table(x$cluster_id[.], x$sample_id[.])
    names(which(rowSums(ns > 30) > 2))
  })
  #ks <- Reduce(intersect, ks)
  ks <- intersect(ks[,1], ks[,2])
  dim(x <- x[, x$cluster_id %in% ks])
  for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])
  table(x$cluster_id, x$sample_id)
  
  # store gene/cell identifiers
  rowData(x)$gene_id <- rownames(x)
  colData(x)$cell_id <- colnames(x)

  
  return(x)
}

