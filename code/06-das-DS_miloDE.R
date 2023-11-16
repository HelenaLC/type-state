suppressPackageStartupMessages({
    library(miloDE)
    library(SingleCellExperiment)
})

fun <- \(x){
    m <- assign_neighbourhoods(x, reducedDim_name = "PCA")
    res <- de_test_neighbourhoods(m, sample_id = "sample_id",
        design = ~ group_id, covariates = c("group_id"))
    res <- na.omit(res)
    if (!is.null(res)) {
        idx <- match(c("gene", "pval", "pval_corrected_across_genes"), names(res))
        names(res)[idx] <- c("gene_id", "p_val", "p_adj")
    }

    return(res)
}