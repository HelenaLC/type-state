#args <- list(list.files("outs", "sco-sim-", full.names=TRUE), "plts/sco-ng.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, select, t, s, b, sim, sco, sco_val)

df <- do.call(res, what=rbind)
th <- seq(0, 1, 0.1)
fd <- lapply(th, \(.) {
    df |>
        mutate(th=.) |>
        filter(!is.na(sco_val)) |>
        group_by(t, s, b, th, sim, sco) |>
        mutate_at("sco_val", \(.) { . <- .-min(.); ./max(.) }) |>
        summarise(n=mean(sco_val <= th), .groups="drop")
}) |> do.call(what=rbind)

gg <- ggplot(fd, 
    aes(th, n, col=sco)) +
    geom_point(size=0.4) +
    facet_grid(t ~ s, labeller=label_both) +
    geom_path(linewidth=0.2, show.legend=FALSE) +
    scale_color_brewer(palette="Paired") +
    guides(col=guide_legend(nrow=2, 
        override.aes=list(size=2, alpha=1))) +
    labs(
        x="0-1 scaled score threshold", 
        y="proportion of features selected") +
    scale_x_continuous(limits=c(0, 1), n.breaks=2) +
    scale_y_continuous(limits=c(0, 1), n.breaks=2) +
    theme_linedraw(9) + theme(
        aspect.ratio=1,
        legend.position="bottom",
        panel.grid.minor=element_blank(),
        legend.key.size=unit(0.5, "lines"),
        strip.text=element_text(color="black"),
        strip.background=element_rect(fill="white"))

ggsave(args[[2]], gg, units = "cm", width = 15, height = 16)
