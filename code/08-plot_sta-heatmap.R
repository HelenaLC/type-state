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
    summarise_at("sta_val", mean, drop = "none")
    #mutate(sta_val = sta_val/max(sta_val, na.rm = TRUE))

gg <- ggplot(df) +
    facet_grid(sel ~ sta) +
    geom_tile(col = "white", aes(t, s, fill = sta_val)) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    labs(x = "type effect", y = "state effect") +
    scale_fill_distiller(NULL,
        palette = "RdBu", na.value = "lightgrey",
        limits = c(-1, 1), n.breaks = 3, direction = 1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA))

ggsave(args[[2]], gg, units = "cm", width = 30, height = 24)
