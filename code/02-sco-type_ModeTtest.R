suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(MKmisc)
})

fun <- \(x) {
    # one vs rest t test
    y <- assay(x, "logcounts")
    ids <- unique(x$cluster_hi)
    res <- sapply(ids, \(k) {
        tmp <- x
        id <- tmp$cluster_id
        j <- !(i <- id == k)
        ij <- c(which(i), which(j))
        df <- data.frame(row.names = ij)
        df[i, "group"] <- "Group1"
        df[j, "group"] <- "Group2"
        df$group <- factor(df$group)
        z <- y[, as.numeric(rownames(df)), drop = FALSE]
        mod.t.test(as.matrix(z), group = df[, "group"])$adj.p.value
    })
    rownames(res) <- rownames(x)
    res <- apply(res, 1, FUN = mean, na.rm = TRUE)

    -log(res)

}
