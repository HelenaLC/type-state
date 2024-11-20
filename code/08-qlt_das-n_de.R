#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-n_de.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(ComplexUpset)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
fdr <- 0.01
df <- lapply(res, select,
    sel, das, gene_id, p_adj) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(sel, das, gene_id) |>
    summarize_at("p_adj", ~any(.x < fdr)) |>
    summarise(n_dge=sum(p_adj), .groups="drop") |>
    mutate(das=gsub("^DS_", "", das), sel=factor(sel, SEL)) 

# plotting
gg <- 
  ggplot(df, 
    aes(das, n_dge, fill=das)) + facet_grid(~sel) +  
    scale_fill_manual(values=c("slateblue", "seagreen", "orange")) +
  ggplot(df, 
    aes(sel, n_dge, fill=sel)) + facet_grid(~das) + 
    scale_fill_manual(values=c("gold", "tomato", "magenta", "cyan", "seagreen1", "seagreen3")) +
plot_layout(nrow=1, 
    widths=c(length(SEL), length(DAS))) &
    labs(x=NULL, y="# DEGs (any adj. p-value < 0.05)") &
    scale_y_continuous(limits=c(0, 2e3), expand=c(0, 0)) &
    plot_annotation(tag_levels="a") &
    geom_col(width=0.75) &
    theme_bw(6) & theme(
        aspect.ratio=2,
        plot.margin=margin(),
        legend.position="none",
        panel.grid.major=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=15, height=4, units="cm")
