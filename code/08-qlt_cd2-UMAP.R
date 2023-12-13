#args <- list(list.files("data/dat/02-rep", "-cd\\.rds$", full.names=TRUE), "plts/dat/cd2-UMAP.pdf")

# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(patchwork)
    library(RColorBrewer)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- do.call(rbind, res)

# wrangling
df <- res[, -grep("PC[0-9]+$", names(res))]
names(df) <- gsub("\\.1", "", names(df))
df$sel <- factor(df$sel, SEL)
df <- df[sample(nrow(df)), ]
if(any(c("X1", "X2") %in% names(df)))
    names(df)[match(c("X1","X2"),names(df))] <- c("UMAP1", "UMAP2")

# aesthetics
aes <- list(
    facet_wrap(~sel, ncol=2, scales="free"),
    geom_point_rast(aes(UMAP1, UMAP2), 
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
gg <- 
    ggplot(df, aes(col=group_id)) +
    scale_color_manual(values=c("royalblue", "tomato")) +
    ggplot(df, aes(col=cluster_id)) +
    scale_color_manual(values=brewer.pal(12, "Paired")) +
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom",
        legend.justification=c(0.5, 1),
        legend.box.spacing=unit(0, "pt"),
        plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=15, height=13, units="cm")
