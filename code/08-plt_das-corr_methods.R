# dependencies
suppressPackageStartupMessages({
    library(poolr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
i <- c("t", "s", "sel", "das", "p_adj", "gene_id")
df <- res |>
    lapply(\(.) .[i]) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(t, s, sel, das, gene_id) |>
    summarize_at("p_adj", ~fisher(.x)$p)
fd <- df |>
    pivot_wider(
        names_from="das", values_from="p_adj",
        id_cols=c("t", "s", "sel", "gene_id")) |>
    group_by(t, s, sel) |>
    summarise(
        .groups="drop",
        edgeR_lemur=cor(DS_edgeR, DS_lemur, method="spearman"),
        edgeR_miloDE=cor(DS_edgeR, DS_miloDE, method="spearman"),
        lemur_miloDE=cor(DS_lemur, DS_miloDE, method="spearman")) |>
    pivot_longer(matches("edgeR|lemur|miloDE")) |>
    mutate(name=gsub("_", "\n", name)) |>
    mutate(values=case_when(value < 0 ~ 0, TRUE ~ value))

# aesthetics
aes <- list(
    facet_grid(sel~name),
    geom_tile(col="white", linewidth=0.1),
    scale_fill_gradientn(
        expression("corr."~X^2), labels=c("<= 0", 1),
        limits=c(0, 1), n.breaks=2, na.value="lightgrey", 
        colors=c("ivory", "pink", "red", "firebrick", "black")),
    scale_x_continuous("type effect", breaks=c(0, 1)),
    scale_y_continuous("state effect", breaks=c(0, 1)),
    coord_equal(expand=FALSE), 
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(vjust=1),
        legend.title=element_text(vjust=1),
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        strip.text.y=element_text(angle=-90),
        plot.tag=element_text(size=9, face="bold")))

# plotting
. <- c(
    "DE", "DEnotDS", "DEgtDS",
    "DS", "DSnotDE", "DSgtDE")
j <- !(i <- fd$sel %in% .)
gg <- 
    ggplot(fd[i, ], aes(t, s, fill=value)) + 
    ggplot(fd[j, ], aes(t, s, fill=value)) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=12.5, units="cm")
