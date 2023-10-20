# args <- list(
#     "code/03-sel-FEAST.R",
#     list.files("outs", "sco-sim-t80,s40,b0,", full.names=TRUE))

source(args[[1]])
sco <- lapply(args[[2]], readRDS)
names(sco) <- sapply(sco, \(.) .$sco[1])

gene_id <- lapply(sco, \(.) .$gene_id)
gene_id <- unique(unlist(gene_id))

df <- data.frame(
    row.names=NULL, wcs, gene_id, 
    sel_val=gene_id %in% fun(sco))

if (!is.null(wcs$sim)) {
    md <- c("t", "s", "b")
    md <- sco$random[1, md]
    df <- data.frame(df, md)
}

saveRDS(df, args[[3]])
