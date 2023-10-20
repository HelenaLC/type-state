x <- c(
    "dplyr",
    "tidyr",
    "edgeR",
    "limma",
    "scran",
    "scater",
    "scuttle",
    "igraph",
    "harmony",
    "splatter",
    "matrixStats",
    "SingleCellExperiment",
    # sco
    "FEAST",
    "scmap",
    "prabhakarlab/DUBStepR",
    # das
    "lemur",
    "miloR",
    "miloDE",
    # sta
    "caret",
    "MLVSBM",
    "cluster",
    "bluster",
    "PCAtools",
    "CellMixS",
    "variancePartition",
    "immunogenomics/lisi",
    # viz
    "UpSetR",
    "ggplot2",
    "ggrastr",
    "patchwork")

# install dependencies
if (!require(BiocManager))
    install.packages("BiocManager")
for (. in x)
    if (!require(., character.only = TRUE))
        BiocManager::install(., ask = FALSE, update = TRUE)

# capture session
for (. in x) {
    . <- gsub(".*/", "", .)
    suppressPackageStartupMessages(
        library(., character.only = TRUE))
}
si <- capture.output(sessionInfo())
writeLines(si, args[[1]])