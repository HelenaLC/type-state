suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(stringr)
    library(ggplot2)
    library(dplyr)
})

# simulated data
res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) %>%
        select(X1, X2, PC1, PC2, sel, 
            t, s, b, group_id, cluster_id, cluster_re)})
df <- do.call(rbind, res)

p1 <- lapply(unique(df$sel), \(x) {
    tmp <- df[df$sel == x, ]
    p <- ggplot(tmp, aes(PC1, PC2, 
        color = factor(cluster_re), shape = factor(group_id))) +
        geom_point(alpha = 0.2, size = 0.6) + 
        facet_grid(t ~ s, labeller = \(.) label_both(.)) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        coord_fixed() + theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank()) + 
        ggtitle(x)
})

# real data
dat <- lapply(args[[1]], \(x) { if (str_detect(x, "dat-")) {
    y <- readRDS(x)
    y %>% select(cluster_re, group_id, PC1, PC2, sel, dat)} 
})

fd <- do.call(rbind, dat)

p2 <- lapply(unique(fd$dat), \(x) {
    tmp <- fd[fd$dat == x, ]
    p <- ggplot(tmp, aes(PC1, PC2, 
        color = factor(cluster_re), shape = factor(group_id))) +
        geom_point(alpha = 0.2, size = 0.6) + 
        facet_wrap(~ sel, ncol = 3, labeller = \(.) label_both(.)) + 
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        coord_fixed() + theme_bw(9) + theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 15),
            panel.grid = element_blank()) + 
        ggtitle(x)
})

plt <- append(p1, p2)


pdf(args[[2]], width = 10, height = 10, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
dev.off()
