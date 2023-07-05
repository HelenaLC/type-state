suppressPackageStartupMessages({
    library(lemur)
    library(SingleCellExperiment)
})

fun <- \(x) {
    fit <- lemur(x, design = ~ group_id)
    fit <- align_harmony(fit) 
    fit <- test_de(fit, 
        contrast = cond(group_id = "Condition1") - cond(group_id = "Condition2"))
    res <- find_de_neighborhoods(fit, group_by = vars(sample_id, group_id))
    res <- res[res$selection == TRUE,]
    idx <- match(c("pval", "adj_pval"), names(res))
    names(res)[idx] <- c("p_val", "p_adj")
    return(res)
}