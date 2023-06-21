#args <- list(list.files("outs", "^sco-.*", full.names = TRUE), "plts/sco.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    mutate(sco_val = case_when(
        sco != "entropy" ~ log10(sco_val),
        .default = sco_val))

# DE genes only (according to 'Splatter')
de <- df[grep("GroupDE", names(df))]
fd <- df[!rowAlls(as.matrix(de) == 1), ]

mg <- sapply(seq_len(ncol(de)), \(i){
    not_i <- setdiff(seq_len(ncol(de)), i)
    mgk <- sapply(not_i, \(j) {
        log(de[,i]/de[,j], base = 2)
    })
    rowMeans(mgk)
})

rownames(mg) <- rownames(de)
# use abs because of cares about down-regulated markers
true <- apply(abs(mg), 1, max)
idx <- which(true > 1)
fc <- df[idx, ] 

gg <- list(
    scale_x_continuous(n.breaks = 3),
    scale_y_continuous(n.breaks = 3),
    geom_density(key_glyph = "point"),
    facet_wrap(~ sco, scales = "free", nrow = 1),
    guides(color = guide_legend(override.aes = list(size = 2))))

p1 <- ggplot(df, aes(sco_val, col = factor(t))) + gg +
    scale_color_brewer(palette = "Blues", "type\neffect")
p2 <- ggplot(fc, aes(sco_val, col = factor(s))) + gg +
    scale_color_brewer(palette = "Reds", "state\neffect")

p3 <- ggplot(df, aes(sco_val, col = factor(t))) + gg +
   scale_color_brewer(palette = "Blues", "type\neffect")
p4 <- ggplot(fc, aes(sco_val, col = factor(s))) + gg +
   scale_color_brewer(palette = "Reds", "state\neffect")

thm <- theme_linedraw(9) + theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    panel.spacing = unit(2, unit = "mm"),
    strip.text = element_text(color = "black", face = "bold"),
    strip.background = element_rect(color = NA, fill = "white"))

plt <- 
    wrap_elements(p1  + plot_layout(guides = "collect") & thm) / 
    wrap_elements(p2  + plot_layout(guides = "collect") & thm) + 
    plot_annotation(tag_levels = "a") &
    theme(
        plot.margin = margin(0, unit = "mm"),
        plot.tag = element_text(size = 9, face = "bold"))

ggsave(args[[2]], plt, units = "cm", width = 30, height = 10)
