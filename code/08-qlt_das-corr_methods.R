#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-corr_methods.pdf")

# dependencies
suppressPackageStartupMessages({
    library(poolr)
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- lapply(res, select,
    sel, das, gene_id, p_adj, dat) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(sel, das, gene_id, dat) |>
    summarize_at("p_adj", ~fisher(.x)$p) |>
    mutate(sel=factor(sel, SEL))

fd <- df |>
    pivot_wider(
        names_from="das", 
        values_from="p_adj",
        id_cols=c("sel", "gene_id", "dat")) |>
    group_by(sel, dat) |>
    summarise(
        .groups="drop",
        edgeR_lemur=cor(DS_edgeR, DS_lemur, 
          method="spearman", use="pairwise.complete.obs"),
        edgeR_miloDE=cor(DS_edgeR, DS_miloDE, 
          method="spearman", use="pairwise.complete.obs"),
        lemur_miloDE=cor(DS_lemur, DS_miloDE, 
          method="spearman", use="pairwise.complete.obs")) |>
    pivot_longer(matches("edgeR|lemur|miloDE")) |>
    mutate(name=gsub("_", "\n", name))

# aesthetics
rng <- range(fd$value, na.rm=TRUE)
rng <- c(
    floor(rng[1]*10)/10, 
    ceiling(rng[2]*10)/10)

# plotting
gg <- lapply(split(fd, fd$dat), \(d) {
  rng <- range(d$value, na.rm=TRUE)
  rng <- c(
    floor(rng[1]*10)/10, 
    ceiling(rng[2]*10)/10)
  ggplot(d, aes(sel, name, fill=value)) + 
    geom_tile(col="white", linewidth=0.1) +
    scale_fill_gradientn(
        expression("Cor("*X[i]^2*","~X[j]^2*")"), 
        limits=rng, breaks=rng, na.value="lightgrey", 
        colors=c("ivory", "pink", "red", "firebrick", "black")) +
    labs(x="feature selection", y="DS methods") +
    coord_equal(2, expand=FALSE) +
    theme_minimal(6) + theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        legend.title=element_text(vjust=1),
        legend.key.size=unit(0.5, "lines"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ggtitle(d$dat[1])}) |> wrap_plots(ncol=3)


# saving
ggsave(args[[2]], gg, width=15, height=4, units="cm")

