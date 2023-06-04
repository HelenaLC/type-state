suppressPackageStartupMessages(
    {
        library(scater)
        library(ggplot2)
    }
)

sce <- readRDS(args[[1]])

name <- basename(args[[1]])
name <- substr(name,4,nchar(name)-4)

if (!is.null(sce)) {
    pca <- plotPCA(sce, 
                   colour_by = "cluster_id", 
                   shape_by = "condition") +
        labs(title = name) + 
        theme(plot.title = element_text(hjust = 0.5)) 
    
    ggsave(args[[2]], plot = pca)
}