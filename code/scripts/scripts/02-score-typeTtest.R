suppressPackageStartupMessages(
    {
        #source("scripts/utils.R")
        library(SummarizedExperiment)
    }
)

fun <- \(x,  
         cluster_to_use = "cluster_id", 
         assay_to_use = "logcounts", 
         penalty = FALSE, 
         fun = mean, 
         fun_pb = "mean"){
    
    # one vs rest t test
    clusters <- colData(x)[,cluster_to_use]
    exprs_mm <- assay(x, assay_to_use)
    n_genes <- 1:nrow(x)
    ttest_per_cluster <- sapply(unique(clusters), \(j){
        temp <- x
        idx <- which(colData(temp)[,cluster_to_use] == j)
        colData(temp)[idx, cluster_to_use] <- 1
        colData(temp)[-idx, cluster_to_use] <- 2
        ttest_res <- sapply(n_genes, function(i){
            if(length(idx) == 1){
                t.test(exprs_mm[i,] ~ 1)$p.value
            }else if(length(idx) > 1 & length(idx) < ncol(temp)){
                t.test(exprs_mm[i,] ~ colData(temp)[,cluster_to_use])$p.value
            }else{
                return(1)
            }
        })
    })
    rownames(ttest_per_cluster) <- rownames(x)
    t_score <- apply(ttest_per_cluster, 1, FUN = fun, na.rm=TRUE)
    
    if(penalty == T){
        colData(x)[,cluster_to_use] <- cluster_ids(x, clustering_to_use)
        pb <- muscat::aggregateData(x, by = cluster_to_use, assay = assay_to_use, fun = fun_pb)
        minmax_norm <- \(x) {
            (x - min(x)) / (max(x) - min(x))
        }
        norm_mat <- apply(assay(pb), 2, minmax_norm)
        max_mat <- apply(norm_mat, 1, max)
        score <- max_mat*(-log(t_score))
        
    }else{
        score <- -log(t_score)
    }
    return(score)
}
