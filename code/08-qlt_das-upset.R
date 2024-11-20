#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-upset.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(UpSetR)
    library(ggplot2)
    library(patchwork)
    library(ComplexUpset)
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
    mutate(das=gsub("^DS_", "", das), sel=factor(sel, SEL))

# split by feature selection method
df_by_sel <- split(df, df$sel)
p1 <- lapply(SEL, \(sel) {
    fd <- df_by_sel[[sel]]
    fd_by_das <- split(fd, fd$das)
    gs <- lapply(fd_by_das, \(.) .$gene_id[.$p_adj])
    # plotting
    thm <- list(legend.text.align=NULL)
    upset(fromList(gs), rev(DAS),
        sort_sets=FALSE, name="",
        height_ratio=1, width_ratio=1/4,
        sort_intersections_by="degree",
        sort_intersections="ascending",
        themes=upset_default_themes(text=element_text(size=6)),
        set_sizes=upset_set_size(position="right") +
            theme(thm, axis.ticks=element_line(linewidth=0.2)) +
            scale_y_continuous(NULL, limits=c(0, 3e3), breaks=c(0, 3e3)),
        base_annotations=list(
            "Intersection size"=
                intersection_size(counts=FALSE) +
                scale_y_continuous(NULL, limits=c(0, 1e3), n.breaks=3) +
                ggtitle(sel) + theme(thm,
                    legend.position="none",
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_line(linewidth=0.2),
                    plot.title=element_text(hjust=0.5, face="bold"))),
        matrix=intersection_matrix(
            geom=geom_point(size=0.8),
            segment=geom_segment(linewidth=0.4)) + theme(thm))
}) |> wrap_plots(nrow=2) & theme(plot.margin=margin(r=5, unit="pt"))

# split by DS analysis method
df_by_das <- split(df, df$das)
p2 <- lapply(DAS, \(das) {
    fd <- df_by_das[[das]]
    fd_by_sel <- split(fd, fd$sel)
    gs <- lapply(fd_by_sel, \(.) .$gene_id[.$p_adj])
    # plotting
    upset(fromList(gs), rev(SEL), 
        sort_sets=FALSE, name="", 
        height_ratio=1, width_ratio=1/4,
        sort_intersections_by="degree",
        sort_intersections="ascending",
        min_size=10, n_intersections=14,
        themes=upset_default_themes(text=element_text(size=6)),
        set_sizes=upset_set_size(position="right") + 
            theme(axis.ticks=element_line(linewidth=0.2)) +
            scale_y_continuous(NULL, limits=c(0, 3e3), breaks=c(0, 3e3)),
        base_annotations=list(
            "Intersection size"=
                intersection_size(counts=FALSE) +
                scale_y_continuous(NULL, limits=c(0, 1200), n.breaks=3) +
                ggtitle(das) + theme(
                    legend.position="none", 
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_line(linewidth=0.2),
                    plot.title=element_text(hjust=0.5, face="bold"))),
        matrix=intersection_matrix(
            geom=geom_point(size=0.4), 
            segment=geom_segment(linewidth=0.2)))
}) |> wrap_plots(nrow=1) & theme(plot.margin=margin(r=5, unit="pt"))

# saving
gg <- 
    wrap_elements(panel=p1) / 
    wrap_elements(panel=p2) + 
    plot_layout(heights=c(2, 1)) +
    plot_annotation(tag_levels="a") &
    theme(
        plot.margin=margin(),
        plot.tag=element_text(size=9, face="bold"))
ggsave(args[[2]], gg, width=15, height=10, units="cm")
