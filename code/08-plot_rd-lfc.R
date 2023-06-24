suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})


res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) 

df <- df[, !duplicated(colnames(df))]
de <- df %>% 
    select(contains("DE"), t, s, b)

gde <- de[grep("GroupDE", names(de))]
mg <- sapply(seq_len(ncol(gde)), \(i){
    not_i <- setdiff(seq_len(ncol(gde)), i)
    mgk <- sapply(not_i, \(j) {
        log(gde[,i]/gde[,j], base = 2)
    })
    rowMeans(mgk)
})

de$mg <- apply(abs(mg), 1, max)


cde <- de[grep("ConditionDE", names(de))]
cg <- sapply(seq_len(ncol(cde)), \(i){
    not_i <- setdiff(seq_len(ncol(cde)), i)
    cgk <- sapply(not_i, \(j) {
        log(cde[,i]/cde[,j], base = 2)
    })
    rowMeans(cgk)
})
de$cg <- apply(cg, 1, max)

gg <- list(
    scale_x_continuous(n.breaks = 3),
    scale_y_continuous(n.breaks = 3),
    geom_density(key_glyph = "point"),
    guides(color = guide_legend(override.aes = list(size = 2))))

p1 <- ggplot(de, aes(x=mg, y = ..scaled.., col = factor(t))) + gg +
    scale_color_brewer(palette = "Blues", "type\neffect") +
    facet_grid( ~ s, labeller = \(.) label_both(.), scales = "free") +
    xlab("GroupDE: max logFC in three groups")

p2 <- ggplot(de, aes(x=cg, y = ..scaled.., col = factor(s))) + gg +
    scale_color_brewer(palette = "Reds", "type\neffect") +
    facet_grid( ~ t, labeller = \(.) label_both(.), scales = "free") +
    xlab("ConditionDE logFC ")

thm <- theme_linedraw(9) + theme(
    panel.grid = element_blank(),
    #axis.title.x = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    panel.spacing = unit(2, unit = "mm"),
    strip.text = element_text(color = "black", face = "bold"),
    strip.background = element_rect(color = NA, fill = "white"))

plt <- wrap_elements(p1  + plot_layout(guides = "collect") & thm) / 
    wrap_elements(p2  + plot_layout(guides = "collect") & thm) + 
    plot_annotation(tag_levels = "a") &
    theme(
        plot.margin = margin(0, unit = "mm"),
        plot.tag = element_text(size = 9, face = "bold"))



ggsave(args[[2]], plt, units = "cm", width = 25, height = 10)


