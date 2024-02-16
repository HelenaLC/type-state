suppressPackageStartupMessages(library(SingleCellExperiment))
fun <- \(x) setNames(rowData(x)$ratio, rownames(x))