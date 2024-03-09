suppressPackageStartupMessages({
    library(lemur)
})

fun <- \(x) {
    ids <- levels(x$group_id)
    fit <- lemur(x, design=~group_id, verbose=FALSE)
    fit <- align_harmony(fit, verbose=FALSE)
    fit <- test_de(fit, cond(group_id=ids[1])-cond(group_id=ids[2]))
    res <- find_de_neighborhoods(fit, vars(sample_id, group_id), 
        verbose=FALSE)
    old <- c("name", "neighborhood", "pval", "adj_pval")
    new <- c("gene_id", "cell_id", "p_val", "p_adj")
    names(res)[match(old, names(res))] <- new
    res$prop <- sapply(res$cell_id, 
      \(i) max(prop.table(table(colData(x)[i, "group_id"]))))
    res$nCells <- sapply(res$cell_id, length)
    cs <- res$cell_id
    cs <- lapply(cs, \(.) paste(.,collapse=";"))
    res$nSubpop <- length(unique(unlist(cs)))
    res$cell_id <- NULL
    return(res)
}