#args <- list(list.files("data/dat/02-rep", "-cd\\.rds$", full.names=TRUE), "plts/dat/cd2-UMAP.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(patchwork)
    library(RColorBrewer)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- lapply(res, select, dat, sel, 
    cluster_id, group_id, UMAP1, UMAP2)

# wrangling
df <- do.call(rbind, res) |>
  mutate(sel=factor(sel, SEL))
df <- df[sample(nrow(df)), ]

# aesthetics
aes <- list(
    facet_wrap(~sel, ncol=2, scales="free"),
    geom_point_rast(aes(UMAP1, UMAP2), 
        shape=16, alpha=0.1, size=0.1),
    guides(col=guide_legend(
        ncol=4, title.position="top",
        override.aes=list(alpha=1, size=2))),
    theme_minimal(6), theme(
        aspect.ratio=1,
        plot.margin=margin(),
        axis.text=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.25, "lines")))

# plotting
ps <- lapply(split(df, df$dat), \(fd) {
  nk <- length(unique(fd$cluster_id))
  ggplot(fd, aes(col=group_id)) +
    scale_color_manual(values=c("royalblue", "tomato")) +
    ggplot(fd, aes(col=cluster_id)) +
    scale_color_manual(values=brewer.pal(nk, "Paired")) +
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
      plot.margin=margin(0),
      legend.position="bottom",
      legend.justification=c(0.5, 1),
      legend.box.spacing=unit(0, "pt"),
      plot.tag=element_text(size=9, face="bold")) 
})

# saving
pdf(args[[2]], onefile=TRUE, width=12/2.54, height=10/2.54)
for (p in ps) print(p); dev.off()
