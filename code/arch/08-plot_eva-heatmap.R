#args <- list.files("outs", "^eva-.*", full.names = TRUE)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
})

res <- lapply(args[[1]], \(x) readRDS(x))
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res) 

lst <- split(df, list(df$sel, df$t, df$s, df$b))
scr <- lapply(lst, \(x) {
    cor_g <- cor(x$cms_g, x$har_g, method = "spearman")
    cor_k <- cor(x$cms_k, x$har_k, method = "spearman")
    dg <- data.frame(sel = x$sel[1], t = x$t[1], 
        s = x$s[1], b = x$b[1], cor = cor_g, col_id = "condition")
    dk <- data.frame(sel = x$sel[1], t = x$t[1], 
        s = x$s[1], b = x$b[1], cor = cor_k, col_id = "cluster")
    rbind(dg, dk)
}) 
res <- do.call(rbind, scr)

gg <- ggplot(res) +
    facet_grid(col_id ~ sel) +
    geom_tile(col = "white", aes(t, s, fill = cor)) +
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

ggsave(args[[2]], gg, units = "cm", width = 30, height = 25)