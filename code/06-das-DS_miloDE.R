suppressPackageStartupMessages({
    library(miloDE)
    library(SingleCellExperiment)
})

fun <- \(x){
    m <- assign_neighbourhoods(x, reducedDim_name = "PCA")
    res <- de_test_neighbourhoods(m, sample_id = "sample_id", 
        design = ~ group_id, covariates = c("group_id"))
    res <- na.omit(res)
    idx <- match(c("pval"), names(res))
    names(res)[idx] <- c("p_val")
    return(res)
}