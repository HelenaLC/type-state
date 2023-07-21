suppressPackageStartupMessages({
    library(UpSetR)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(grid)
})


sel <- lapply(args[[1]], \(x) { if (!str_detect(x,"dat")) readRDS(x) })
names(sel) <- sapply(sel, \(.) .$sel[1])
df <- do.call(rbind, sel) 
df$sim_gene <- paste0(df$sim, df$gene_id)


# sim data
lst <- split(df, df$sel)
sel_idx <- lapply(lst, \(x) x$sim_gene[x$sel_val == TRUE])
p1 <- upset(fromList(sel_idx),
    order.by = "freq",
    nsets = length(sel_idx),
    set_size.show = TRUE,
    scale.sets = "identity",
    mainbar.y.label = "Intersection Size - Simulated data",
    text.scale = 2)

# real data
les <- lapply(args[[1]], \(x) { if (str_detect(x,"dat")) readRDS(x) })
names(les) <- sapply(sel, \(.) .$sel[1])
fd <- do.call(rbind, les)

plt <- lapply(unique(fd$dat), \(x){
    df <- fd[fd$dat == x,]
    lst <- split(df, df$sel)
    sel_idx <- lapply(lst, \(x) x$gene_id[x$sel_val == TRUE])
    p <- upset(fromList(sel_idx),
        order.by = "freq",
        nsets = length(sel_idx),
        set_size.show = TRUE,
        scale.sets = "identity",
        mainbar.y.label = paste0("Intersection Size - ", x),
        text.scale = 2)
    return(p)
})
names(plt) <- unique(fd$dat)
plt[["Simulated data - all"]] <- p1

pdf(args[[2]], width = 20, height = 8, onefile = TRUE)
for (i in names(plt)) {
    print(plt[[i]])
}
dev.off()

#ggsave(filename = args[[2]], plot = p, width = 10, height = 8)