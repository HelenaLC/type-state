#args <- list(list.files("outs", "sco-", full.names=TRUE), "plts/sco-ng.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, select, t, s, b, sim, sco, sco_val)

df <- do.call(res, rbind)
th <- seq(0, 1, 0.1)
fd <- lapply(th, \(.) {
    df |>
        mutate(th=.) |>
        filter(!is.na(sco_val)) |>
        group_by(t, s, b, th, sim, sco) |>
        mutate_at("sco_val", \(.) { . <- .-min(.); ./max(.) }) |>
        summarise(n=mean(sco_val <= th), .groups="drop")
}) |> do.call(what=rbind)

lab <- labeller(.rows=label_value, .cols=label_both)

p1 <- ggplot(fd, aes(th, n, col=t, group=t)) + 
    facet_grid(sco~s, labeller=lab) +
    scale_color_distiller("type\neffect", palette="Blues", n.breaks=2)

p2 <- ggplot(fd, aes(th, n, col=s, group=s)) + 
    facet_grid(sco~t, labeller=lab) +
    scale_color_distiller("state\neffect", palette="Reds", n.breaks=2)

gg <- wrap_plots(p1, p2, nrow=1) +
    plot_annotation(tag_levels="a") &
    geom_path(linewidth=0.4) &
    labs(
        x="0-1 scaled score threshold", 
        y="proportion of features selected") &
    scale_x_continuous(limits=c(0, 1), n.breaks=2) &
    scale_y_continuous(limits=c(0, 1), n.breaks=2) &
    theme_linedraw(6) + theme(
        aspect.ratio=1,
        legend.position="bottom",
        panel.grid.minor=element_blank(),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        plot.tag=element_text(face="bold"),
        legend.title.align=1,
        legend.title=element_text(vjust=1.5),
        strip.text=element_text(color="black"),
        strip.background=element_rect(fill="white"))

ggsave(args[[2]], gg, units = "cm", width = 18, height = 13)
