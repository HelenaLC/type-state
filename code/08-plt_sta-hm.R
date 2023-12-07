args <- list(list.files("outs/sim", "^sta-.*", full.names=TRUE), "plts/sim/sta-hm.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- do.call(rbind, res)
.f <- \(df) df %>%
    group_by(t, s, b, sel, sta) %>%
    summarise_at("sta_val", mean) %>%
    mutate(sta_val=case_when(
        sta_val < 0 ~ 0, 
        .default=sta_val))

# separate ground truth selections
gt <- c(
    "DE", "DEnotDS", "DEgtDS",
    "DS", "DSnotDE", "DSgtDE")
df_des <- .f(filter(df,  sel %in% gt))
df_sel <- .f(filter(df, !sel %in% gt))

# aesthetics
aes <- list(
    facet_grid(sel ~ sta),
    geom_tile(
        aes(t, s, fill=sta_val), 
        col="white", linewidth=0.1),   
    scale_fill_gradientn(NULL,
        colors=c("ivory", "gold", "red", "navy"),
        na.value="lightgrey", limits=c(0, 1), n.breaks=2),
    labs(x="type effect", y="state effect"),
    scale_x_continuous(breaks=c(0, 1)),
    scale_y_continuous(breaks=c(0, 1)),
    coord_fixed(expand=FALSE),
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(vjust=1),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        plot.tag=element_text(size=9, face="bold")))

# plotting
gg <- ggplot(df_des) + ggplot(df_sel) +
    plot_layout(ncol=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=18, units="cm")
