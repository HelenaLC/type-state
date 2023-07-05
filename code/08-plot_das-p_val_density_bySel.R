# args <- list(
#     list.files("outs", "^das.*", full.names = TRUE),
#     "plts/das-p_val_density.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(stringr)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(.) .[c("sel","t", "s", "das", "p_val", "p_adj")])

typ <- vapply(res, \(.) .$das[1], character(1))
ds <- !(da <- grepl("^DA", typ))
da <- do.call(rbind, res[da])
ds <- do.call(rbind, res[ds])
all <- do.call(rbind, res)
all <- all %>% 
    mutate(
        DAS = ifelse(str_detect(das, "DS"), "DS", "DA"),
        clustering = ifelse(str_detect(das, "limma|edgeR"), "Yes", "No"))

plt <- lapply(unique(all$sel), \(s){
    df <- all[all$sel == s,]
    gg <- ggplot(df, aes(p_val, after_stat(ndensity), col = das, linetype = DAS)) +
        facet_grid(s ~ t, labeller = \(.) label_both(.)) +
        geom_density(key_glyph = "point") +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        scale_x_continuous("p-value", breaks = seq(0, 1, 0.5), expand = c(0, 0.05)) +
        scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
        ggtitle(s) + 
        theme_bw(9) + theme(
            panel.grid = element_blank(),
            #legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.key.size = unit(0.5, "lines"),
            axis.text.x = element_text(angle = 45, hjust = 1))
})

pdf(args[[2]], width = 10, height = 8, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
dev.off()



