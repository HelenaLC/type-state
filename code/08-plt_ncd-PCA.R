#args <- list(list.files("data/02-rep", "-cd", full.names=TRUE), "plts/ncd-PCA.pdf")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(stringr)
    library(ggplot2)
    library(dplyr)
})

# simulated data
idx <- grepl("^sim", basename(args[[1]]))
res <- lapply(args[[1]][idx], readRDS)
res <- lapply(res, select, 
    PC1, PC2, sel, t, s, b, 
    group_id, cluster_id, cluster_re)
df <- do.call(rbind, res)

p1 <- lapply(split(df, df$sel), \(fd) 
    ggplot(fd, aes(PC1, PC2, 
        color=factor(cluster_re), shape=factor(group_id))) +
        geom_point(alpha=0.2, size=0.6) + 
        facet_grid(t ~ s, labeller=\(.) label_both(.)) +
        scale_x_continuous(n.breaks=3) +
        scale_y_continuous(n.breaks=3) +
        ggtitle(fd$sel[1]) +
        coord_fixed() + 
        theme_bw(9) + theme(
            legend.position="none",
            panel.grid=element_blank()))

# real data
idx <- grepl("^dat", basename(args[[1]]))
res <- lapply(args[[1]][idx], readRDS)
res <- lapply(res, select, 
    PC1, PC2, sel, dat, 
    group_id, cluster_re)
df <- do.call(rbind, res)

p2 <- lapply(split(df, df$dat), \(fd)
    ggplot(fd, aes(PC1, PC2, 
        color=factor(cluster_re), 
        shape=factor(group_id))) +
        geom_point(alpha=0.2, size=0.6) + 
        facet_wrap(~ sel, ncol=3, labeller=\(.) label_both(.)) + 
        scale_x_continuous(n.breaks=3) +
        scale_y_continuous(n.breaks=3) +
        ggtitle(fd$dat[1]) + 
        coord_fixed() + 
        theme_bw(9) + theme(
            legend.position="none",
            plot.title=element_text(hjust=0.5, size=15),
            panel.grid=element_blank()))

pdf(args[[2]], width=10, height=10, onefile=TRUE)
print(p1); print(p2); dev.off()
