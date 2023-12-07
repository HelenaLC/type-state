#args <- list(list.files("data/sim/01-pro", "-rd\\.rds", full.names=TRUE), "plts/sim/rd1-lfc.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- do.call(rbind, res) 
de <- grep("GroupDE", names(df))
ds <- grep("ConditionDE", names(df))

df$de <- sapply(de, \(i) {
    sapply(setdiff(de, i), \(j) 
        log2(df[[i]]/df[, j])) |>
        rowMeans() }) |> 
    abs() |> rowMeans()

df$ds <- sapply(ds, \(i) {
    sapply(setdiff(ds, i), \(j) 
        log2(df[[i]]/df[, j])) |>
        rowMeans() }) |> 
    abs() |> rowMeans()

fd <- df |>
    pivot_longer(all_of(c("de", "ds"))) |>
    mutate_at("name", factor, c("de", "ds"), c("between\ngroups", "between\nconditions"))

# aesthetics
aes <- list(
    geom_density(linewidth=0.4, key_glyph="point"),
    scale_x_continuous("mean |logFC|", n.breaks=3),
    scale_y_continuous("scaled density", n.breaks=2),
    guides(color=guide_legend(override.aes=list(size=1))),
    theme_bw(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        panel.spacing=unit(2, unit="mm"),
        legend.key.size=unit(0.25, "lines"),
        strip.text=element_text(color="black"),
        strip.background=element_rect(color=NA, fill="white")))
.labs <- labeller(name=label_value, s=label_both, t=label_both)
fac_s <- facet_grid(name ~ s, labeller=.labs, scales="free")
fac_t <- facet_grid(name ~ t, labeller=.labs, scales="free")
pal_s <- scale_color_brewer(palette="Reds", "state\neffect", limits=seq(0, 1, 0.2))
pal_t <- scale_color_brewer(palette="Blues", "type\neffect", limits=seq(0, 1, 0.2))

# plotting
p1 <- ggplot(fd, aes(value, ..scaled.., col=factor(t))) + pal_t + fac_s + aes
p2 <- ggplot(fd, aes(value, ..scaled.., col=factor(s))) + pal_s + fac_t + aes
gg <- (p1 + p2) +
    plot_layout(ncol=1) +
    plot_annotation(tag_levels="a") &
    theme(
        plot.margin=margin(),
        plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=15, height=8, units="cm")
