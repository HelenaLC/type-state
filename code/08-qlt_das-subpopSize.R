#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-nCluster.pdf")

# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggrastr)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- lapply(res, select, 
    das, dat, sel, prop, nCells) |>
    do.call(what=rbind) |>
    distinct(das, dat, sel, nCells, prop) |>
    mutate(sel=factor(sel, SEL))

gg <- ggplot(df, aes(das, nCells, col=das)) + 
  geom_violin_rast(position = position_dodge(width = 0.8), trim = FALSE) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position = "none") +
  facet_grid2(sel~dat, scales="free", independent = "y") +
  scale_color_brewer(palette = "Paired")

ggsave(args[[2]], gg, width=8, height=20, units="cm")