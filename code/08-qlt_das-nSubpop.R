#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-nCluster.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- lapply(res, select, 
  das, dat, sel, nSubpop) |>
  do.call(what=rbind) |>
  distinct(das, dat, sel, nSubpop)


gg <- ggplot(df, aes(sel, das, fill=as.numeric(nSubpop))) + 
    geom_tile(col="white") +
    facet_grid(~dat, scales="free", labeller=\(.) label_both(.)) +
    labs(x="feature selection", y="DS methods")  + 
    geom_text(aes(label = nSubpop), size = 2.5) +
    theme_bw() +
    theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        legend.title=element_text(vjust=1),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        scale_fill_distiller(NULL,
        palette="YlOrRd", direction = 1)


ggsave(args[[2]], gg, width=20, height=8, units="cm")