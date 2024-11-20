#args <- list(list.files("outs/dat", "^das-.*Kang", full.names=TRUE), "plts/dat/das-bar.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidytext)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
fdr <- 0.05
df <- lapply(res, select,
    sel, das, gene_id, p_adj) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(sel, das, gene_id) |>
    summarize_at("p_adj", ~any(.x < fdr)) |>
    filter(p_adj) |> select(-p_adj) |>
    mutate(
        sel=factor(sel, SEL),
        das=factor(gsub("^DS_", "", das), DAS))

# plotting
ex <- expansion(0, 0.6)
ey <- expansion(0.02, 0)

fd <- df |>
    group_by(sel, gene_id) |>
    count() |>
    filter(n == length(DAS)) |>
    group_by(sel) |>
    count()

pal_das <- scale_fill_manual(values=c("slateblue", "seagreen", "orange"))
pal_sel <- scale_fill_manual(values=c("gold", "tomato", "magenta", "cyan", "seagreen1", "seagreen3"))

p1 <- ggplot(fd, aes(reorder(sel, n), n, fill=sel)) + 
    labs(x="selection method") +
    ggtitle("all DS methods") +
    scale_y_continuous(limits=c(0, 1200), expand=ey) +
    pal_sel

fd <- df |>
    group_by(das, gene_id) |>
    count() |>
    filter(n == length(SEL)) |>
    group_by(das) |>
    count()

p2 <- ggplot(fd, aes(reorder(das, n), n, fill=das)) + 
    labs(x="DS analysis method") +
    ggtitle("all selections") +
    scale_y_continuous(NULL, limits=c(0, 1200), expand=ey) +
    pal_das

fd <- df |>
    filter(sel == "random") |>
    group_by(das) |>
    count() 

p3 <- ggplot(fd, aes(reorder(das, n), n, fill=das)) + 
    labs(x="DS analysis method") +
    ggtitle("random only") +
    scale_y_continuous(NULL, limits=c(0, 2e3), expand=ey) +
    pal_das

ws <- c(length(SEL), length(DAS), length(DAS))
gg <- wrap_plots(p1, p2, p3, widths=ws) +
    plot_layout(guides="collect") &
    labs(y=sprintf("# DEGs (any adj. p-value < %s)", fdr)) &
    geom_col(width=0.8) &
    scale_x_reordered(expand=ex) &
    theme_minimal(6) & theme(
        panel.grid=element_blank(), 
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.5, "lines"),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[2]], gg, width=10, height=6, units="cm")
