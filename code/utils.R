.ei <- \(x, 
         sample = "sample_id",
         group = "condition")
{
    ids <- unique(x[[sample]])
    idx <- match(ids, x[[sample]])
    ei <- colData(x)[idx, c(sample, group)]
    ei <- data.frame(ei, row.names = NULL)
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

.ds <- \(x, 
  cluster = names(cluster_codes(x))[1],
  sample = "sample_id", group = "condition")
{
  # setup model & contrast matrix
  ei <- .ei(x, sample, group)
  dm <- createDesignMatrix(ei, group)
  cm <- createContrast(c(0, 1))
  rowData(x)$marker_class <- "state"
  # run differential abundance (DA) analysis
  se <- diffcyt(x, ei, 
    clustering_to_use = cluster, 
    design = dm, contrast = cm, verbose = FALSE, 
    analysis_type = "DS", method_DS = "diffcyt-DS-limma")$res
  data.frame(rowData(se))
}



.all_score <- \(x, dim){
  x <- cluster(x, 
               xdim = dim, ydim = dim,
               features = rownames(x), 
               seed = seed, verbose = FALSE)
  se <- .se(x,  "cluster_id", "sample_id", "exprs", "median")
  se$condition <- x$condition[match(se$sample_id, x$sample_id)]
  da <- .da(x)
  ds <- .ds(x)
  res <- .score(se, da, ds)
  res$marker_class <- marker_classes(x)[res$marker_id]
  return(res)
}




.calculate_stability <- \(x, clustering_to_use = names(cluster_codes(sce)[1]), score, weighted = FALSE){
  if(weighted == T){
    assay(x, "exprs") <- score$type_score*assay(x, "exprs") 
  }
  clusters <- cluster_codes(sce)[clustering_to_use][,1]
  stability <- lapply(clusters, \(j){
    idx <- which(cluster_ids(x, clustering_to_use)==j)
    temp_x <- x[,idx]
    exprs <- assay(temp_x, "exprs")
    exprs[exprs < 0] <- 0
    apply(exprs, 1, function(x) prop.table(table(exprs)))
    w_j <- rowSums(exprs==0)/ncol(exprs)
    sigma_j <- apply(data.frame(exprs), 1, var, na.rm=TRUE)
    mean_j <- apply(data.frame(exprs), 1, mean, na.rm=TRUE)
    n_genes <- nrow(exprs)
    s_j <- 1-(w_j+sigma_j/mean_j)/n_genes
    #data.frame(w = w_j, sigma = sigma_j, mean = mean_j)
  })
  res <- matrix(unlist(stability), ncol = length(stability))
  rownames(res) <- names(stability[[1]])
  return(res)
}

.roc_type <- \(x,
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
