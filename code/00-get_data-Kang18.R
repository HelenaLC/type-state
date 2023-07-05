suppressPackageStartupMessages({
    library(ExperimentHub)
    library(SingleCellExperiment)
    library(muscat)
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
    
    x$id <- paste0(x$stim, x$ind)
    (x <- prepSCE(x, 
        kid = "cell", # subpopulation assignments
        gid = "stim",  # group IDs (ctrl/stim)
        sid = "id",   # sample IDs (ctrl/stim.1234)
        drop = TRUE))
    
    colData(x) <- DataFrame(
        cluster_id = x$cell,
        sample_id = paste0(x$stim, x$ind),
        group_id = x$stim)
    for (. in names(colData(x)))
        x[[.]] <- factor(x[[.]])
    
    # store gene/cell identifiers
    rowData(x)$gene_id <- rownames(x)
    colData(x)$cell_id <- colnames(x)
    
}