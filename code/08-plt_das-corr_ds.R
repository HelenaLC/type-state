# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- lapply(res, select, t, s, 
    sel, das, p_adj, gene_id) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(t, s, sel, das, gene_id) |>
    summarize_at("p_adj", min)

das <- paste(DAS, collapse="|")
fd <- df |>
    pivot_wider(
        names_from="das", values_from="p_adj",
        id_cols=c("t", "s", "sel", "gene_id")) |>
    rename_with(\(.) gsub("^DS_", "", .), matches(das)) |>
    group_by(t, s, sel) |>
    summarise(
        .groups="drop",
        muscat_lemur=cor(muscat, lemur, method="spearman"),
        muscat_miloDE=cor(muscat, miloDE, method="spearman"),
        lemur_miloDE=cor(lemur, miloDE, method="spearman")) |>
    pivot_longer(matches(das)) |>
    mutate(
        name=gsub("_", "\n", name),
        value=case_when(value < 0 ~ 0, TRUE ~ value))

fd <- mutate(fd, sel=factor(sel, c(DES, SEL)))
j <- !(i <- fd$sel %in% DES)
df_des <- fd[i, ]
df_sel <- fd[j, ]

# aesthetics
aes <- list(
    facet_grid(sel~name),
    scale_fill_gradientn(
        "corr. min(adj. p-value)",
        colors=c("ivory", "pink", "red", "firebrick", "black"),
        na.value="lightgrey", labels=c("<= 0", 1), limits=c(0, 1), n.breaks=2),
    geom_tile(col="white", linewidth=0.1, aes(t, s, fill=value)),
    scale_x_continuous("type effect", breaks=c(0, 1)),
    scale_y_continuous("state effect", breaks=c(0, 1)),
    coord_equal(expand=FALSE), 
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        legend.title=element_text(vjust=1),
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        strip.text.y=element_text(angle=-90),
        plot.tag=element_text(size=9, face="bold")))

# plotting
gg <- ggplot(df_des) + ggplot(df_sel) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=12.5, units="cm")
