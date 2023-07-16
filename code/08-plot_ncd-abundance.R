suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(stringr)
    library(ggplot2)
    library(dplyr)
})

# simulated data
res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) })
df <- do.call(rbind, res)

plt <- lapply(unique(df$sel), \(x) {
    tmp <- df[df$sel == x, ]
    p <- ggplot(tmp, aes(x = sample_id, fill = factor(cluster_re))) +
        geom_bar(position = "fill") + 
        facet_grid(t ~ s, labeller = \(.) label_both(.), scales = "free") +
        scale_y_continuous(n.breaks = 3) +
        theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) + 
        ggtitle(x)
})

pdf(args[[2]], width = 10, height = 10, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
dev.off()