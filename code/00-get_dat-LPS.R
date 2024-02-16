suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
})

fun <- \() {
  # load SCE
  x <- readRDS("data/SCE_annotation.rds")
  
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
  for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])
  
  # store gene/cell identifiers
  rowData(x)$gene_id <- rownames(x)
  colData(x)$cell_id <- colnames(x)
  
  return(x)
}

