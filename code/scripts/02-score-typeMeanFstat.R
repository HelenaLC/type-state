suppressPackageStartupMessages(
    {
        library(FEAST)
        library(SummarizedExperiment)
    }
)

fun <- \(x, 
         assay_to_use = "logcounts", 
         cluster_to_use = "cluster_id",
         sample_to_use = "sample_id"){
    
    F_stat_per_sample <- sapply(unique(colData(x)[,sample_to_use]), \(i){
        temp <- x[,which(colData(x)[,sample_to_use] == i)]
        Y <- as.matrix(assay(temp, assay_to_use))
        F_stat <- cal_F2(Y, colData(temp)[,cluster_to_use])
        F_stat$F_scores
    })
    rownames(F_stat_per_sample) <- rownames(x)
    rowMeans(F_stat_per_sample)
}