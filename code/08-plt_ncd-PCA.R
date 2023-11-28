#args <- list(list.files("data/02-rep", "-cd", full.names=TRUE), "plts/ncd-PCA.pdf")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(ggplot2)
    library(dplyr)
    library(ggrastr)
})

# simulated data
idx <- grepl("^sim", basename(args[[1]]))
res <- lapply(args[[1]][idx], readRDS)
res <- lapply(res, select, 
    PC1, PC2, sel, t, s, b, 
    group_id, cluster_id, cluster_re)
df <- do.call(rbind, res)

# symmetrically round axis limits to nearest 5
.f <- \(.) round(max(abs(.))/5)*5
pm <- c(-1, 1)
xs <- pm*.f(df$PC1)
ys <- pm*.f(df$PC2)

# setup more easily interpretable color palette
df$id <- paste(df$cluster_id, df$group_id)
n <- 7; i <- c(2, 4)
pal <- c(
    hcl.colors(n, "reds")[i], 
    hcl.colors(n, "blues")[i], 
    hcl.colors(n, "greens")[i])
names(pal) <- levels(factor(df$id))

p1 <- lapply(split(df, df$sel), \(fd) 
    ggplot(
        fd[sample(nrow(fd)), ], 
        aes(PC1, PC2, col=id)) +
        geom_point_rast(shape=16, alpha=0.2, size=0.1) + 
        facet_grid(t ~ s, labeller=\(.) label_both(.)) +
        guides(col=guide_legend(nrow=2, 
            override.aes=list(alpha=1, size=2))) +
        ggtitle(fd$sel[1]) + 
        scale_color_manual(NULL, values=pal) +
        scale_x_continuous(
            "1st principal component", 
            n.breaks=3, limits=xs) +
        scale_y_continuous(
            "2nd principal component", 
            n.breaks=3, limits=ys) +
        theme_bw(6) + theme(
            legend.position="bottom",
            panel.grid.minor=element_blank(),
            axis.title=element_text(hjust=0),
            legend.key.size=unit(0.5, "lines")))

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

pdf(args[[2]], width=8, height=8, onefile=TRUE)
print(p1); print(p2); dev.off()
