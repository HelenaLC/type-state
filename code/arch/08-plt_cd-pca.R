#args <- list(list.files("data/01-fil", "cd\\.rd*", full.names=TRUE), "plts/cd-pca.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
})

res <- lapply(args[[1]], readRDS)
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

gg <- ggplot(
    df[sample(nrow(df)), ], 
    aes(PC1, PC2, col=id)) +
    geom_point_rast(shape=16, alpha=0.2, size=0.1) + 
    facet_grid(t ~ s, labeller=\(.) label_both(.)) +
    guides(col=guide_legend(nrow=2, 
        override.aes=list(alpha=1, size=2))) +
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
        legend.key.size=unit(0.5, "lines"))

ggsave(args[[2]], gg, units="cm", width=15, height=15)
