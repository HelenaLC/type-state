suppressPackageStartupMessages({
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
    }
    return(res)
}