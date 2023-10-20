#args <- list(list.files("data/02-rep", "-cd", full.names=TRUE), "plts/ncd-abundance.pdf")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(stringr)
    library(ggplot2)
    library(dplyr)
})

idx <- grepl("^sim-", basename(args[[1]]))
res <- lapply(args[[1]][idx], readRDS)
res <- lapply(res, select, t, s, sel, cluster_id, cluster_re)
df <- do.call(rbind, res)

p1 <- lapply(split(df, df$sel), \(fd)
    ggplot(fd, aes(x=cluster_id, 
        fill=factor(cluster_re))) +
        geom_bar(position="fill") + 
        facet_grid(t ~ s, labeller=\(.) label_both(.), scales="free") +
        scale_y_continuous(n.breaks=3) +
        ggtitle(fd$sel[1]) +
        theme_bw(9) + theme(
            legend.position="none",
            panel.grid=element_blank(),
            axis.text.x=element_text(angle=45, vjust=1, hjust=1)))

idx <- grepl("^dat-", basename(args[[1]]))
res <- lapply(args[[1]][idx], readRDS)
res <- lapply(res, select, dat, sel, cluster_id, cluster_re)
df <- do.call(rbind, res)

p2 <- lapply(split(df, df$dat), \(fd)
    ggplot(fd, aes(x=cluster_id, 
        fill=factor(cluster_re))) +
        geom_bar(position="fill") + 
        facet_wrap(~ sel, ncol=3, labeller=\(.) label_both(.)) + 
        scale_y_continuous(n.breaks=3) +
        ggtitle(fd$dat[1]) +
        theme_bw(9) + theme(
            legend.position="none",
            panel.grid=element_blank(),
            axis.text.x=element_text(angle=45, vjust=1, hjust=1)))

pdf(args[[2]], width=10, height=10, onefile=TRUE)
print(p1); print(p2); dev.off()