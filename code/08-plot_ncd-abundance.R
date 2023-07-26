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
    p <- ggplot(tmp, aes(x = cluster_id, fill = factor(cluster_re))) +
        geom_bar(position = "fill") + 
        facet_grid(t ~ s, labeller = \(.) label_both(.), scales = "free") +
        scale_y_continuous(n.breaks = 3) +
        theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
        ggtitle(x)
})

out <- lapply(args[[1]], \(x) { if (str_detect(x, "dat-")) {
    y <- readRDS(x)
    y %>% select(cluster_re, cluster_id, sel, dat)} 
})
fd <- do.call(rbind, out)

p2 <- lapply(unique(fd$dat), \(x) {
    tmp <- fd[fd$dat == x, ]
    p <- ggplot(tmp, aes(x = cluster_id, fill = factor(cluster_re))) +
        geom_bar(position = "fill") + 
        facet_wrap(~ sel, ncol = 3, labeller = \(.) label_both(.)) + 
        scale_y_continuous(n.breaks = 3) +
        theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
        ggtitle(x)
})

plt <- append(p1, p2)

pdf(args[[2]], width = 10, height = 10, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
dev.off()