suppressPackageStartupMessages({
    library(Seurat)
    library(DUBStepR)
    library(SingleCellExperiment)
})

fun <- \(x) {
    y <- GetAssayData(ScaleData(as.Seurat(x)))
    z <- DUBStepR(y, optimise.features=TRUE)
    sco_val <- z$corr.info$corr.range
    names(sco_val) <- z$corr.info$feature.genes
    sel_val <- names(sco_val) %in% z$optimal.feature.genes
    data.frame(
        row.names=NULL, sco_val, sel_val,
        gene_id=gsub("-", "_", names(sco_val)))
}