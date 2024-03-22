#args <- list(list.files("data/dat/01-pro", "\\.rds$", full.names=TRUE), "plts/dat/cd2-UMAP.pdf")

# dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
  library(RColorBrewer)
  library(ggpubr)
  library(SingleCellExperiment)
})

# loading
res <- lapply(args[[1]], \(x) {
  sce <- readRDS(x)
  data.frame(reducedDim(sce, "HARMONY")[, c("PC1", "PC2")], 
    cluster_id=sce$cluster_id, group_id=sce$group_id, dat = basename(x)
    )
})
df <- do.call(rbind, res)
df$dat <- gsub("\\.rds$", "", df$dat)

# aesthetics
aes <- list(
  facet_wrap(~dat, ncol=2, scales="free"),
  geom_point_rast(aes(PC1, PC2), 
    shape=16, alpha=0.1, size=0.1),
  guides(col=guide_legend(
    ncol=4, title.position="top",
    override.aes=list(alpha=1, size=1))),
  theme_minimal(6), theme(
    aspect.ratio=1,
    plot.margin=margin(),
    axis.text=element_blank(),
    panel.grid=element_blank(),
    axis.title=element_text(hjust=0),
    panel.border=element_rect(fill=NA),
    legend.key.size=unit(0.25, "lines")))

# plotting
gg <- lapply(split(df, df$dat), \(fd) {
  nCells <- length(unique(fd$cluster_id))
  ggplot(fd, aes(col=group_id)) +
    scale_color_manual(values=c("royalblue", "tomato")) +
    ggplot(fd, aes(col=cluster_id)) +
    scale_color_manual(values=get_palette("Paired",k=nCells)) +
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(fd$dat[1],tag_levels="a") &
    aes & theme(
      plot.margin=margin(),
      legend.position="bottom",
      legend.justification=c(0.5, 1),
      legend.box.spacing=unit(0, "pt"),
      plot.tag=element_text(size=7, face="bold")) 
  
}) |> wrap_plots(ncol=1)


# saving
ggsave(args[[2]], gg, width=15, height=6, units="cm")