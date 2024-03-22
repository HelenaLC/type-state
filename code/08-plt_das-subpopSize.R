#args <- list(list.files("outs/sim", "^das-", full.names=TRUE), "plts/dat/das-nCluster.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- lapply(res, select, 
    das, t, s, sel, prop, nCells) |>
    do.call(what=rbind) |>
    distinct(das, t, s, sel, nCells, prop) |>
    group_by(das, t, s, sel) |>
    summarise_at("nCells", median) |>
    mutate(nCells=case_when(nCells > 1500 ~ 1500, TRUE ~ nCells)) 



aes <- list(
    geom_tile_rast(
        aes(t, s, fill=nCells), 
        col="white", linewidth=0.1),
    facet_grid(sel~das),
    scale_fill_gradientn(
        "median # of cells/subpopulation",
        colors=c("ivory", "gold", "red", "navy"), 
        na.value="lightgrey", limits=range(df$nCells),
        breaks=c(95, 800, 1500),
        labels=c("95", "800", ">= 1500")),
    labs(x="type effect", y="state effect"),
    scale_x_continuous(breaks=c(0, 1)),
    scale_y_continuous(breaks=c(0, 1)),
    coord_fixed(expand=FALSE),
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.title=element_text(vjust=1.5),
        plot.tag=element_text(size=9, face="bold")))

df <- mutate(df, sel=factor(sel, c(DES, SEL)))
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# plotting
gg <- ggplot(df_des) + ggplot(df_sel) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=12.5, units="cm")