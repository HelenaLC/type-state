suppressPackageStartupMessages(
    {
        library(scater)
        library(ggplot2)
    }
)

sce <- readRDS(args[[1]])

name <- basename(args[[1]])
name <- substr(name,1,nchar(name)-8)

if (!is.null(sce)) {
    umap <- plotUMAP(sce, colour_by = "cluster_id", shape_by = "condition") +
        labs(title = name) + 
        theme(plot.title = element_text(hjust = 0.5)) 
    
    ggsave(args[[2]], plot = umap)
}


