#args <- list(list.files("outs", "^sta-.*", full.names = TRUE), "plts/sta.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    group_by(t, s, b, sel, sta) %>%
    summarise_at("sta_val", mean)

gg <- ggplot(df) +
    facet_grid(sta ~ sel) +
    geom_tile(col = "white", aes(t, s, fill = sta_val)) +
    scale_x_continuous(breaks = seq(0, 1, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    labs(x = "type effect", y = "state effect") +
    scale_fill_gradientn("statistic", 
        limits = c(-1, 1), n.breaks = 3, 
        colors = hcl.colors(11, "Spectral")) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(args[[2]], gg, units = "cm", width = 15, height = 12)
