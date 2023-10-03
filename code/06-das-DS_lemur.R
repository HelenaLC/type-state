suppressPackageStartupMessages({
    library(lemur)
    library(SingleCellExperiment)
})

fun <- \(x) {
    fit <- lemur(x, design = ~ group_id)

    #fit <- align_harmony(fit) # got error, why?
    #fit <- align_by_grouping(fit, "cluster_id")

    cd1 <- unique(x$group_id)[1]
    cd2 <- unique(x$group_id)[2]
    fit <- test_de(fit, 
        contrast = cond(group_id = cd1) - cond(group_id = cd2))
    res <- find_de_neighborhoods(fit, group_by = vars(sample_id, group_id))
    res <- res[res$selection == TRUE,]
    idx <- match(c("pval", "adj_pval"), names(res))
    names(res)[idx] <- c("p_val", "p_adj")
    return(res)
}