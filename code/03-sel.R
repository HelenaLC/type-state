source(args[[1]])
sco <- lapply(args[[2]], readRDS)

sco <- lapply(sco, \(.) {
    if (length(unique(.$sco)) > 1) 
        split(., .$sco) else list(.)
})
sco <- unlist(sco, recursive=FALSE)
names(sco) <- sapply(sco, \(.) .$sco[1])

gene_id <- lapply(sco, \(.) .$gene_id)
gene_id <- unique(unlist(gene_id))

df <- data.frame(
    row.names=NULL, wcs, gene_id, 
    sel_val=gene_id %in% fun(sco))

if (!is.null(wcs$sim)) {
    md <- c("t", "s", "b")
    md <- sco[[1]][1, md]
    df <- data.frame(df, md)
}

saveRDS(df, args[[3]])
