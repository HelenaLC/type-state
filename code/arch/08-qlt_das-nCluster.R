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

df <- lapply(res, \(fd) {
    if ("cluster_re" %in% colnames(fd)) 
        nSubpopulation <- max(fd$cluster_re)
    else if ("Nhood" %in% colnames(fd))
        nSubpopulation <- max(fd$Nhood)
    else 
        nSubpopulation <- length(unique(fd[,"cell_id"]))
    data.frame(fd[,c("das", "dat", "sel")], nSubpopulation=nSubpopulation)
})

df <- do.call(rbind, df) |> distinct(das, dat, sel, nSubpopulation)

gg <- ggplot(df, aes(sel, das, fill=as.numeric(nSubpopulation))) + 
    geom_tile(col="white") +
    facet_grid(~dat, scales="free") +
    labs(x="feature selection", y="DS methods")  + 
    geom_text(aes(label = nSubpopulation), size = 3) +
    theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        legend.title=element_text(vjust=1),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey")

ggsave(args[[2]], gg, width=20, height=8, units="cm")