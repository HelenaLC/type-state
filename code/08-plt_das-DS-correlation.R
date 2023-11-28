# args <- list(
#     list.files("outs", "^das.*", full.names=TRUE),
#     "plts/das-p_val_density.pdf")

suppressPackageStartupMessages({
    library(poolr)
    library(dplyr)
    library(tidyr)
    library(ggpubr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
})

# DS results
idx <- grep("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
i <- c("t", "s", "sel", "das", "p_adj", "gene_id")
df <- res |>
    lapply(\(df) df[i]) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(t, s, sel, das, gene_id) |>
    summarize_at("p_adj", ~-log(fisher(.x)$p))

# ground truth
rd <- list.files("data/01-fil", "*rd.rds", full.names=TRUE)
rd <- do.call(rbind, lapply(rd, readRDS))

# wrangling
de <- grep("^GroupDE", names(rd))
ds <- grep("^ConditionDE", names(rd))

# type effect
rd$de <- sapply(de, \(i) {
    j <- setdiff(de, i)
    lfc <- log2(rd[[i]]/rd[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

# state effect
rd$ds <- sapply(ds, \(i) {
    j <- setdiff(ds, i)
    lfc <- log2(rd[[i]]/rd[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

# merge with DS results
rd <- rd[c("t", "s", "gene_id", "de", "ds")]
df <- merge(x=df, y=rd, by=c("t", "s", "gene_id"))

# correlation between DS results & ground truth
dfs <- split(df, df[. <- c("t", "s", "sel", "das")])
res <- lapply(dfs, \(df) {
    t <- cor(df$p_adj, df$de, method="spearman")
    s <- cor(df$p_adj, df$ds, method="spearman")
    data.frame(
        row.names=NULL, df[1, .], 
        cor=c("cor_t", "cor_s"), 
        cor_val=c(t, s))
}) |> 
    do.call(what=rbind) |>
    na.omit()

lab <- "Spearman's correlation with %s effect"
p1 <- ggplot(res[res$cor == "cor_t", ]) + ggtitle(sprintf(lab, "type"))
p2 <- ggplot(res[res$cor == "cor_s", ]) + ggtitle(sprintf(lab, "state"))

ps <- lapply(list(p1, p2), \(p) p + facet_grid(das ~ sel) +
        geom_tile(col="white", aes(t, s, fill=cor_val)) +
        scale_x_continuous("type effect", n.breaks=2) +
        scale_y_continuous("state effect", n.breaks=2) +
        scale_fill_distiller(NULL,
            palette="RdYlBu", na.value="lightgrey",
            limits=c(-1, 1), n.breaks=3, direction=-1) +
        coord_fixed(expand=FALSE) + 
        theme_minimal(9) + theme(
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_rect(fill=NA),
            legend.key.width=unit(1, "lines"),
            legend.key.height=unit(0.5, "lines")))

# wrangling
fd <- df |>
    # group_by(t, s, sel, das) |>
    # mutate(
    #     p_adj=p_adj-min(p_adj),
    #     p_adj=p_adj/max(p_adj)) |>
    # ungroup() |>
    mutate(
        across(all_of(c("t", "s")), 
            ~factor(., sort(unique(.))))) |>
    pivot_wider(
        names_from="das", values_from="p_adj",
        id_cols=setdiff(names(df), c("das", "p_adj"))) |>
    pivot_longer(
        names_to="effect",
        cols=c("de", "ds")) |>
    mutate(effect=factor(effect, 
        levels=c("de", "ds"), 
        labels=c("type", "state"))) |>
    ungroup() |>
    arrange(value)

pa <- ggplot(fd, aes(DS_edgeR, DS_miloDE, col=value)) 
pb <- ggplot(fd, aes(DS_edgeR, DS_lemur, col=value))

p3 <- wrap_plots(pa, pb, ncol=1) +
    plot_layout(guides="collect") & 
    facet_grid(effect~sel) &
    geom_point_rast(shape=16, alpha=0.1, size=0.5) &
    stat_cor(method="pearson", aes(label=..r.label..)) &
    scale_color_viridis_c() &
    scale_x_continuous(n.breaks=3) &
    scale_y_continuous(n.breaks=3) &
    theme_linedraw(6) & theme(
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        legend.key.size=unit(0.5, "lines"))

pdf(args[[2]], width=10, height=7, onefile=TRUE)
for (p in c(ps, list(p3))) print(p); dev.off()