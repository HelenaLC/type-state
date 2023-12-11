#args <- list(list.files("outs/sim", "^sco-", full.names=TRUE), "plts/sim/sco-ecdf.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- res|>
    do.call(what=rbind) |>
    mutate(sco=factor(sco, SCO))
df$sco_t <- paste(df$sco, df$t, sep="_")
df$sco_s <- paste(df$sco, df$s, sep="_")
de <- grep("GroupDE", names(df))
ds <- grep("ConditionDE", names(df))
fd <- df[rowAnys(df[de] != 1) & !rowAnys(df[ds] != 1), ]
df <- bind_rows(.id="sub", 
    list(all=df, DEnotDS=fd)) |>
    mutate(sub=factor(sub, c("all", "DEnotDS")))

# aesthetics
aes <- list(
    scale_x_continuous(n.breaks=3),
    scale_y_continuous(n.breaks=2),
    labs(x="score value", y="ECDF"),
    facet_grid(sub ~ sco, scales="free_x"),
    stat_ecdf(linewidth=0.4, key_glyph="point"),
    guides(color=guide_legend(override.aes=list(size=1))),
    theme_bw(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.spacing=unit(2, unit="mm"),
        legend.key.size=unit(0.25, "lines"),
        strip.text=element_text(color="black"),
        strip.background=element_rect(color=NA, fill="white")))
pal_s <- scale_color_brewer(palette="Reds", "state\neffect", limits=seq(0, 1, 0.2))
pal_t <- scale_color_brewer(palette="Blues", "type\neffect", limits=seq(0, 1, 0.2))

# plotting
p1 <- ggplot(df, aes(sco_val, group=sco_t, col=factor(t))) + pal_t
p2 <- ggplot(df, aes(sco_val, group=sco_s, col=factor(s))) + pal_s
gg <- (p1 + p2)+ plot_layout(ncol=1) &
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=15, height=8, units="cm")
