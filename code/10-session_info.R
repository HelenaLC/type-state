x <- c(
    "dplyr",
    "tidyr",
    "limma",
    "scran",
    "scater",
    "scuttle",
    "igraph",
    "harmony",
    "splatter",
    "matrixStats",
    "SingleCellExperiment",
    # das
    "poolr",
    "edgeR",
    "lemur",
    "miloR",
    "MarioniLab/miloDE",
    # sta
    "caret",
    "MLVSBM",
    "bluster",
    "cluster",
    "CellMixS",
    "variancePartition",
    "immunogenomics/lisi",
    # viz
    "UpSetR",
    "ggplot2",
    "ggrastr",
    "patchwork",
    "ComplexUpset")

# install dependencies
if (!require(BiocManager))
    install.packages("BiocManager")
for (. in x)
    if (!require(basename(.), character.only = TRUE))
        BiocManager::install(., ask = FALSE, update = TRUE)

# capture session
for (. in x) {
    . <- gsub(".*/", "", .)
    suppressPackageStartupMessages(
        library(., character.only = TRUE))
}
si <- capture.output(sessionInfo())
writeLines(si, args[[1]])