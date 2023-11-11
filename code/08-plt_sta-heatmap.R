#args <- list(list.files("outs", "^sta-.*", full.names=TRUE), "plts/sta-heatmap.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    group_by(t, s, b, sel, sta) %>%
    summarise_at("sta_val", mean) %>%
    mutate(sta_val=case_when(
        sta_val < 0 ~ 0, 
        .default=sta_val))

gg <- ggplot(df, 
    aes(t, s, fill=sta_val)) + 
    facet_grid(sta ~ sel) +
    geom_tile(col="white") +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        limits=c(0, 1), n.breaks=3, direction=-1) +
    labs(x="type effect", y="state effect") +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) +
    coord_fixed(expand=FALSE) + 
    theme_minimal(6) + theme(
        legend.position="bottom",
        panel.grid=element_blank(),
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"))

f <- \(.) length(unique(.))
h <- (w <- 15)*f(df$sta)/f(df$sel)
ggsave(args[[2]], gg, units="cm", width=w, height=h)
