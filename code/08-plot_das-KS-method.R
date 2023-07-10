suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x)})
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(x) {x %>% select(t, s, b, das, sel, p_val) })

df <- do.call(rbind, res) %>%
    group_by(t, s, b, das, sel) %>%
    filter(!is.na(p_val)) %>%
    summarise(
        ks_stat = ifelse(n() >= 2, ks.test(p_val, "punif")$statistic, NA),
        p_value = ifelse(n() >= 2, ks.test(p_val, "punif")$p.value, NA))

gg <- ggplot(df) +
    facet_grid(t ~ s, labeller = \(.) label_both(.)) +
    geom_tile(col = "white", aes(sel, das, fill = ks_stat)) +
    #scale_x_continuous(breaks = c(0, 1)) +
    #scale_y_continuous(breaks = c(0, 1)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    labs(x = "feature selection", y = "DA/DS method") +
    scale_fill_distiller(NULL,
        palette = "Blues", na.value = "lightgrey",
        limits = c(0, 1), n.breaks = 3, direction = 1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA))
        #axis.text.x = element_text(angle = 75, vjust = 0.5, hjust = 1))

ggsave(args[[2]], gg, units = "cm", width = 25, height = 18)