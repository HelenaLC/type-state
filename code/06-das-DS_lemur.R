suppressPackageStartupMessages({
    library(lemur)
})

fun <- \(x) {
        ids <- levels(x$group_id)
        fit <- lemur(x, design=~group_id)
        fit <- test_de(fit, contrast=cond(group_id=ids[1])-cond(group_id=ids[2]))
        res <- find_de_neighborhoods(fit,
            group_by=vars(sample_id, group_id),
            useNames = TRUE)
        old <- c("name", "indices", "pval", "adj_pval")
        new <- c("gene_id", "cell_id", "p_val", "p_adj")
        names(res)[match(old, names(res))] <- new
        return(res)
}