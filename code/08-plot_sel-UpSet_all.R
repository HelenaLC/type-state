suppressPackageStartupMessages({
    library(UpSetR)
    library(ggplot2)
    library(dplyr)
})


sel <- lapply(args[[1]], readRDS)
names(sel) <- sapply(sel, \(.) .$sel[1])
df <- do.call(rbind, sel) 

lst <- split(df, df$sel)
sel_idx <- lapply(lst, \(x) x$gene_id[x$sel_val == TRUE])
p <- upset(fromList(sel_idx),
    order.by = "freq",
    nsets = length(sel_idx),
    set_size.show = TRUE,
    scale.sets = "identity")

pdf(args[[2]], width = 10, height = 8)
p
dev.off()

#ggsave(filename = args[[2]], plot = p, width = 10, height = 8)