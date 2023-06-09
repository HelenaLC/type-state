suppressPackageStartupMessages({
    library(scater)
    library(ggplot2)
})

sce <- readRDS(args[[1]])

lab <- basename(args[[1]])
lab <- substr(lab, 4, nchar(lab)-4)

if (!is.null(sce)) {
    plt <- plotPCA(sce, 
        colour_by = "cluster_id", shape_by = "group_id") +
        coord_equal() + ggtitle(lab) + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) 
    ggsave(args[[2]], plot = plt)
}