suppressPackageStartupMessages({
  library(MouseGastrulationData)
  library(SingleCellExperiment)
  library(scater)
  library(scDblFinder)
})

fun <- \() {
  # load SCE
  x <- WTChimeraData()
  
  rownames(x) <- uniquifyFeatureNames(rowData(x)$ENSEMBL, rowData(x)$SYMBOL)
  # remove undetected genes
  x <- x[rowSums(counts(x)) > 0, ]
  
  colData(x) <- DataFrame(
    cluster_id=x$celltype.mapped,
    sample_id=x$sample,
    group_id=x$tomato)
  
  # drop multiplet cells
  x <- x[,x$cluster_id!="Doublet"]
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
