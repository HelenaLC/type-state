#args <- list(list.files("outs", "^sco-.*", full.names = TRUE), "plts/sco.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)

gg <- list(
    geom_density(key_glyph = "point"),
    facet_wrap(~ sco, scales = "free", nrow = 1),
    guides(color = guide_legend(override.aes = list(size = 2))),
    theme_linedraw(9), theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        strip.text = element_text(color = "black", face = "bold"),
        strip.background = element_rect(color = NA, fill = "white")))

p1 <- ggplot(df, aes(sco_val, col = factor(t))) + gg +
    scale_color_brewer(palette = "Blues", "type\neffect")
p2 <- ggplot(df, aes(sco_val, col = factor(s))) + gg +
    scale_color_brewer(palette = "Reds", "state\neffect")

plt <- p1 / p2

ggsave(args[[2]], plt, units = "cm", width = 15, height = 8)
