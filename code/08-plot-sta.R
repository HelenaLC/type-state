#args <- list(list.files("outs", "^sta-.*", full.names = TRUE), "plts/sta.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res)
df$id <- with(df, paste(sco, sel, sep = ","))

fd <- df %>%
    group_by(t, s, b, sco, sel, sta) %>%
    summarise_at("sta_val", mean)

gg <- ggplot(fd) +
    facet_grid(sta ~ sco) +
    geom_tile(col = "white", aes(t, s, fill = sta_val)) +
    scale_x_continuous(breaks = seq(0, 1, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    scale_fill_distiller("statistic", 
        limits = c(0, 1), n.breaks = 2, 
        palette = "Blues", direction = 1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(args[[2]], gg, units = "cm", width = 15, height = 12)
