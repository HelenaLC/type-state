suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(SingleCellExperiment)
  library(variancePartition)
  library(bluster)
})

sco <- lapply(args[[3]], readRDS)
sce <- readRDS(args[[2]])
sel <- unlist(strsplit(wcs$sel, "_"))
num <- round((as.numeric(wcs$num)/10)*nrow(sce))

source(args[[1]])

sco <- lapply(sco, \(.) {
  if (length(unique(.$sco)) > 1) 
    split(., .$sco) else list(.)
})
sco <- unlist(sco, recursive=FALSE)
names(sco) <- sapply(sco, \(.) .$sco[1])

if (length(sel)==2) {
  t <- rank(sco[[sel[1]]]$sco_val)
  s <- rank(sco[[sel[2]]]$sco_val)
  o <- order(t-s, decreasing=TRUE)
} else {
  t <- rank(sco[[sel[1]]]$sco_val)
  o <- order(t, decreasing=TRUE)
}
res <- rowData(sce)$gene_id[o[seq_len(num)]]
rowData(sce)$sel_val <- sel_val <- rowData(sce)$gene_id %in% res
sce <- runPCA(sce, subset_row=sel_val)
sce$cluster_re <- clusterCells(sce,
  use.dimred = "PCA",
  BLUSPARAM = NNGraphParam(cluster.fun = "louvain",
    cluster.args = list(resolution = 0.4)))

sta <- fun(sce)
ex <- names(df <- sta)
wcs <- wcs[setdiff(names(wcs), ex)]
df <- data.frame(wcs, 
  df, 
  nGenes = num)
df <- rbind(df, data.frame(dat=wcs$dat,
  num=wcs$num,
  sel=wcs$sel,
  sta="nCluster",
  sta_val=as.numeric(length(levels(sce$cluster_re))),
  nGenes = num))
saveRDS(df, args[[4]])

