# args <- list(
#     list.files("outs/sim", "^das-", full.names=TRUE),
#     "plts/sim/das-KS_effects.pdf")

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
df <- lapply(res, select, 
    t, s, sel, das, p_val) |>
    do.call(what=rbind) |>
    group_by(t, s, sel, das) |>
    filter(!is.na(p_val)) |>
    summarise(.groups="drop",
        ks_stat=ifelse(n() >= 2, ks.test(p_val, "punif")$statistic, NA),
        p_value=ifelse(n() >= 2, ks.test(p_val, "punif")$p.value, NA)) |>
    mutate(ks_stat=case_when(ks_stat > 0.2 ~ 0.2, TRUE ~ ks_stat)) |>
    mutate(das=gsub("^DS_", "", das))

# split selection
des <- c(
    "DE", "DEnotDS", "DEgtDS",
    "DS", "DSnotDE", "DSgtDE")
j <- !(i <- df$sel %in% des)
fd <- df[j, ]; df <- df[i, ]

# aesthetics
aes <- list(
    geom_tile(
        col="white", linewidth=0.1, 
        aes(t, s, fill=ks_stat)),
    facet_grid(das ~ sel, labeller=\(.) label_value(.)),
    labs(x="type effect", y="state effect"),
    scale_x_continuous(breaks=c(0, 1)),
    scale_y_continuous(breaks=c(0, 1)),
    scale_fill_gradientn(
        "KS stat. p-values\nuniformly distributed", 
        guide=guide_colorbar(reverse=TRUE),
        colors=c("black", "red", "gold", "white"),
        na.value="lightgrey", limits=c(0, 0.2),
        n.breaks=2, labels=c("0", ">= 0.2")),
    coord_fixed(expand=FALSE),
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(vjust=1.5),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines")))

# plotting
gg <- ggplot(df) / ggplot(fd) +
    plot_layout(guides="collect") &
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(), 
        legend.position="bottom",
        plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=12, height=14.5, units="cm")
