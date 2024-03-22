#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-nCluster.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(poolr)
    library(ggh4x)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- lapply(res, select, 
    das, dat, sel, gene_id, p_adj) |>
    do.call(what=rbind) |>
    group_by(das, dat, sel, gene_id) |>
    summarise_at("p_adj", min) |>
    #summarize_at("p_adj", ~fisher(.x)$p) |>
    group_by(das, dat, sel) |>
    summarize(nDE = as.numeric(sum(p_adj < 0.05)), .groups = "drop") |>
    mutate(sel=factor(sel, SEL))


gg <- ggplot(df, aes(sel,nDE,col=sel)) +
  geom_bar(stat="identity", fill="white", position=position_dodge()) +
  facet_grid2(dat~das, scales="free", independent="y") +
  scale_color_brewer(palette = "Paired") +
  theme_minimal(6) & theme(
    legend.position="bottom",
    legend.justification=c(0.5, 1),
    legend.box.spacing=unit(0, "pt"),
    panel.grid.minor=element_blank(), 
    panel.border=element_rect(fill=NA),
    plot.tag=element_text(size=9, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[2]], gg, width=15, height=7, units="cm")