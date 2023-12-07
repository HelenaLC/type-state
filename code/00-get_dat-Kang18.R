suppressPackageStartupMessages({
    library(ExperimentHub)
    library(SingleCellExperiment)
})

fun <- \() {
    # load SCE from EH
    eh <- ExperimentHub()
    q <- query(eh, "Kang18_8vs8")
    x <- eh[[q$ah_id]]
    
    # remove undetected genes
    x <- x[rowSums(counts(x)) > 0, ]
    
    # drop unassigned & multiplet cells
    x <- x[, !is.na(x$cell)]
    x <- x[, x$multiplets == "singlet"]
    
    # break sample pairing
    . <- sort(unique(x$ind))
    j <- setdiff(., i <- .[seq_len(4)])
    y <- x[, 
        x$ind %in% i & x$stim == "ctrl" |
        x$ind %in% j & x$stim == "stim" ]
    ncol(y)/ncol(x)
    table(y$ind, y$stim)
    
    colData(y) <- DataFrame(
        cluster_id=y$cell,
        sample_id=y$ind,
        group_id=y$stim)
    for (. in names(colData(y)))
        y[[.]] <- factor(y[[.]])
    
    # store gene/cell identifiers
    rowData(y)$gene_id <- rownames(y)
    colData(y)$cell_id <- colnames(y)
    
    return(y)
}