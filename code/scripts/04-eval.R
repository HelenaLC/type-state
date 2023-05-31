suppressPackageStartupMessages(
    {
        source("code/scripts/utils.R")
        library(bluster)
        library(CellMixS)
        library(SingleCellExperiment)
        library(dplyr)
        library(Metrics)
    }
)

sce <- readRDS(args[[1]])
output <- args[[2]]
name <- basename(args[[1]])
name <- substr(name,1,nchar(name)-8)
sel <- readRDS(paste0("output/sel/", name, "_sel.rds"))
#sim <- args[[3]]

sil <- .silhouette(sce)
np <- .neighborPurity(sce)
cms <- .cms_cluster(sce, k = 50, fun = mean)
auc <- .auc(sce, sel_idx = sel)
res <- list(silhouette = sil, 
            cms = cms, 
            neighborPurity = np, 
            precision = auc$precision, 
            recall = auc$recall, 
            accuracy = auc$accuracy)

saveRDS(res, output)