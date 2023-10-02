#args <- list(list.files("data/01-fil", "cd\\.rd*", full.names = TRUE), "plts/pca.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

df <- do.call(rbind, lapply(args[[1]], readRDS))
# df <- df %>% rename("PC1" = "1",
#     "PC2" = "2")

df$id <- paste(df$cluster_id, df$group_id)
pal <- c(
    hcl.colors(8, "reds")[c(2, 4)], 
    hcl.colors(8, "blues")[c(2, 4)], 
    hcl.colors(8, "greens")[c(2, 4)])
names(pal) <- levels(factor(df$id))

gg <- ggplot(df, aes(PC1, PC2, col = id)) +
    geom_point(shape = 16, alpha = 0.2, size = 0.2) + 
    facet_grid(t ~ s, labeller = \(.) label_both(.)) +
    scale_color_manual(values = pal) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    coord_fixed() + theme_bw(9) + theme(
        legend.position = "none",
        panel.grid = element_blank())

ggsave(args[[2]], gg, units = "cm", width = 15, height = 14)
