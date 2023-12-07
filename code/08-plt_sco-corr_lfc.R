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
df <- do.call(rbind, res) 

# wrangling
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
    group_by(t, s, sco) |>
    summarise(
        .groups="drop",
        de=cor(sco_val, de, method="spearman"),
        ds=cor(sco_val, ds, method="spearman")) |>
    pivot_longer(all_of(c("de", "ds"))) |>
    mutate(name=factor(name, c("de", "ds"), c("type", "state"))) |>
    mutate(sco=factor(sco, c("random", "HVG", 
        sort(unique(grep("type", df$sco, value=TRUE))),
        sort(unique(grep("state", df$sco, value=TRUE))))))

# aesthetics
rng <- range(fd$value, na.rm=TRUE)
rng <- c(
    floor(rng[1]*10)/10, 
    ceiling(rng[2]*10)/10)

# plotting
gg <- ggplot(fd, aes(t, s, fill=value)) +
    geom_tile(col="white", linewidth=0.1) +
    facet_grid(name~sco) + 
    scale_fill_gradient2(
        "corr. |mean logFC|", limits=rng, breaks=c(rng, 0), 
        na.value="lightgrey", low="blue", mid="grey95", high="red") +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) +
    labs(x="type effect", y="state effect") +
    coord_fixed(expand=FALSE) +
    theme_minimal(6) + theme(
        plot.margin=margin(),
        legend.position="bottom",
        legend.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        legend.title=element_text(vjust=1),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        panel.border=element_rect(fill=NA, linewidth=0.2))

# saving
ggsave(args[[2]], gg, width=11, height=5, units="cm")
