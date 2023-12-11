#args <- list(list.files("data/sim/02-rep", "-cd\\.rds", full.names=TRUE), "plts/sim/cd2-pca.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
})

# loading
res <- lapply(args[[1]], readRDS)

# wrangling
df <- res |>
    do.call(what=rbind) |>
    mutate(
        k=gsub("Group", "cluster", cluster_id),
        g=gsub("Condition", "group", group_id),
        i=paste(k, g), sel=factor(sel, c(DES, SEL)))

# aesthetics
n <- 7; i <- c(2, 4)
pal <- c(
    hcl.colors(n, "reds")[i], 
    hcl.colors(n, "blues")[i], 
    hcl.colors(n, "greens")[i])
names(pal) <- levels(factor(df$i))
.f <- \(.) round(max(abs(.))/5)*5
pm <- c(-1, 1)
xs <- pm*.f(df$PC1)
ys <- pm*.f(df$PC2)

# plotting
ps <- lapply(split(df, df$sel), \(fd) 
    ggplot(fd[sample(nrow(fd)), ], aes(PC1, PC2, col=i)) +
        scale_color_manual(NULL, values=pal) +
        geom_point(shape=16, alpha=0.2, size=0.1) + 
        facet_grid(t ~ s, labeller=\(.) label_both(.)) +
        guides(col=guide_legend(nrow=2, override.aes=list(alpha=1, size=2))) +
        scale_x_continuous("1st principal component", n.breaks=3, limits=xs) +
        scale_y_continuous("2nd principal component", n.breaks=3, limits=ys) +
        ggtitle(fd$sel[1]) + theme_bw(6) + theme(
            legend.position="bottom",
            panel.grid.minor=element_blank(),
            axis.title=element_text(hjust=0),
            legend.key.size=unit(0.5, "lines")))

# saving
pdf(args[[2]], onefile=TRUE, width=15/2.54, height=16/2.54)
for (p in ps) print(p); dev.off()
