suppressPackageStartupMessages(library(SingleCellExperiment))
fun <- \(x) setNames(rowData(x)$bio, rownames(x))