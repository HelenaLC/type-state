x <- c(
    "bluster",
    "CellMixS",
    "SingleCellExperiment"
)

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