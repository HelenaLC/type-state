# args <- list(
#     list.files("outs", "^das.*", full.names = TRUE),
#     "plts/das-p_val_density.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(.) .[c("t", "s", "das", "p_val", "p_adj")])

typ <- vapply(res, \(.) .$das[1], character(1))
ds <- !(da <- grepl("^DA", typ))
da <- do.call(rbind, res[da])
ds <- do.call(rbind, res[ds])

gg <- ggplot(ds, aes(p_val, after_stat(ndensity), col = das)) +
    facet_grid(s ~ t, labeller = \(.) label_both(.)) +
    geom_density(key_glyph = "point") +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.5), expand = c(0, 0.05)) +
    scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
    theme_bw(9) + theme(
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(args[[2]], gg,
    units = "cm", width = 15, height = 9,
    dpi = 300, useDingbats = FALSE)
