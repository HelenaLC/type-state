suppressPackageStartupMessages({
    library(scmap)
})

fun <-  \(x) {
    rowData(x)$feature_symbol <- rownames(x)
    x <- selectFeatures(x, suppress_plot = TRUE)
    rowData(x)$scmap_scores
}
