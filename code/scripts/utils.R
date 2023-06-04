.ei <- \(x, 
         sample = "sample_id",
         group = "condition")
{
    ids <- unique(x[[sample]])
    idx <- match(ids, x[[sample]])
    ei <- colData(x)[idx, c(sample, group)]
    ei <- data.frame(ei, row.names = NULL)
} 


.one_vs_rest_t_test <- \(x,  cluster_to_use = "cluster_id", 
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
        cell1 <- which(colData(temp)[,cluster_to_use] == j)
        cell2 <- which(colData(temp)[,cluster_to_use] != j)
        group.info <- data.frame(row.names = c(cell1, cell2))
        group.info[cell1, "group"] <- "Group1"
        group.info[cell2, "group"] <- "Group2"
        group.info[, "group"] <- factor(x = group.info[, "group"])
        data.use <- exprs_mm[, as.numeric(rownames(group.info)), drop = FALSE]
        ttest_res <- sapply(n_genes, function(i){
            if (length(unique(group.info[, "group"])) == 2) {
                t.test(data.use[i,] ~ group.info[, "group"])$p.value
            } else {
                t.test(data.use[i,] ~ 1)$p.value
            }
        })
    })
    rownames(ttest_per_cluster) <- rownames(x)
    t_score <- apply(ttest_per_cluster, 1, FUN = fun, na.rm=TRUE)
    
    if(penalty == T){
        colData(x)[,cluster_to_use] <- cluster_ids(x, clustering_to_use)
        pb <- muscat::aggregateData(x, 
                                    by = cluster_to_use, 
                                    assay = assay_to_use, 
                                    fun = fun_pb)
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

.one_vs_rest_wilcox_score <- \(x, 
                               assay_to_use = "logcounts", 
                               fun = "mean", 
                               cluster_to_use = "cluster_id", 
                               penalty = FALSE) {
    exprs_mm <- assay(x, assay_to_use)
    n_genes <- 1:nrow(x)
    wilcox_per_cluster <- sapply(unique(colData(x)[, cluster_to_use]), \(j){
        temp <- x
        cell1 <- which(colData(temp)[,cluster_to_use] == j)
        cell2 <- which(colData(temp)[,cluster_to_use] != j)
        group.info <- data.frame(row.names = c(cell1, cell2))
        group.info[cell1, "group"] <- "Group1"
        group.info[cell2, "group"] <- "Group2"
        group.info[, "group"] <- factor(x = group.info[, "group"])
        data.use <- exprs_mm[, as.numeric(rownames(group.info)), drop = FALSE]
        wilcox_res <- sapply(n_genes, function(i){
            if (length(unique(group.info[, "group"])) == 2) {
                wilcox.test(data.use[i,] ~ group.info[, "group"])$p.value
            } else {
                wilcox.test(data.use[i,] ~ 1)$p.value
            }
        })
        
    })
    rownames(wilcox_per_cluster) <- rownames(x)
    mean_wilcox <- apply(wilcox_per_cluster, 1, mean, na.rm=TRUE)
    
    if(penalty){
        pb <- muscat::aggregateData(x, 
                                    by = cluster_to_use, 
                                    assay = assay_to_use, 
                                    fun = fun)
        minmax_norm <- \(x) {
            (x - min(x)) / (max(x) - min(x))
        }
        norm_mat <- apply(assay(pb), 2, minmax_norm)
        max_mat <- apply(norm_mat, 1, max)
        score <- max_mat*(-log(mean_wilcox))
        
    }else{
        score <- -log(mean_wilcox)
    }
    
}

.cms_cluster <- \(x, 
                  sample = "sample_id", 
                  cluster = "cluster_id", 
                  fun = mean, ...){
    
    res <- sapply(unique(colData(x)[, sample]), \(i){
        temp <- x[,which(colData(x)[, sample] == i)]
        dim <- dim(reducedDim(temp, "PCA"))[2]
        cms.res <- cms(temp, group = cluster, n_dim = dim, ...)
        return(fun(cms.res$cms))
    })
    fun(res)
}

.neighborPurity <- \(x, cluster = "cluster_id"){
    pure.res <- neighborPurity(reducedDim(x, "PCA"), colData(x)[, cluster])
    pure <- mean(pure.res$purity)
}


.silhouette <- \(x, 
                 sample_id = "sample_id", 
                 cluster_id = "cluster_id", 
                 fun = mean, 
                 dr = "PCA"){
    sil <- sapply(unique(colData(x)[,sample_id]), \(i){
        idx <- which(colData(x)[,sample_id] == i)
        temp <- x[, idx]
        if(length(unique(colData(temp)[,cluster_id])) > 1){
            sil.approx <- approxSilhouette(reducedDim(temp, dr), clusters=colData(temp)[,cluster_id])
            sil.data <- as.data.frame(sil.approx)
            sil.data$closest <- factor(ifelse(sil.data$width > 0, sil.data$cluster, sil.data$other))
            fun(sil.data$width)
        }else{
            return(NA)
        }
    })
    return(fun(sil))
}

.auc <- \(x,
          sel_idx = sel_idx){
    groupDE <- data.frame(rowData(x)) %>%
        select(GroupDE.Group1, GroupDE.Group2, GroupDE.Group3) 
    true_markers <- groupDE[rowSums(groupDE == 1) == ncol(groupDE) - 1, ]
    true_markers_idx <- match(rownames(true_markers), rownames(x))
    predicted <- true <- rep(0, nrow(x))
    predicted[sel_idx] <- 1
    true[true_markers_idx] <- 1
    recall <- recall(predicted, true)
    precision <- precision(predicted, true)
    accuracy <- accuracy(predicted, true)
    return(list(recall = recall, precision = precision, accuracy = accuracy))
}


.createContrast <- function(contrast) {
    
    contrast_matrix <- matrix(contrast, ncol = 1)
    
    contrast_matrix
}

.createDesignMatrix <- function(experiment_info, cols_design = NULL) {
    
    stopifnot(any(class(experiment_info) %in% c("data.frame", "tbl_df", "tbl")) || is(experiment_info, "DataFrame"))
    experiment_info <- as.data.frame(experiment_info)
    
    # terms for design matrix
    if (is.character(cols_design)) {
        stopifnot(all(cols_design %in% colnames(experiment_info)))
        terms <- cols_design
    } else if (is.numeric(cols_design) | is.logical(cols_design)) {
        terms <- colnames(experiment_info)[cols_design]
    } else if (is.null(cols_design)) {
        # default: all columns
        terms <- colnames(experiment_info)
    }
    
    # create design matrix
    formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
    
    design <- model.matrix(formula, data = experiment_info)
    
    design
}

