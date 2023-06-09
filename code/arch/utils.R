

.se <- \(x, 
  # variables to group by (entropy 
  # is computed across i for every j)
  i = "cluster_id", 
  j = "sample_id",    
  assay = "exprs",  # measurement data to use
  fun = "median",   # summary function used for aggregation
  base = 2) 
{
  # compute pseudo-bulks by cluster-sample
  y <- aggregateAcrossCells(x, 
    colData(x)[c(i, j)], coldata.merge = FALSE,
    use.assay.type = assay, statistics = fun)
  # split cell indices by samples
  idx <- split(seq(ncol(y)), y[[j]])
  # compute proportions across clusters
  # (matrix of dim. features x clusters)
  p <- lapply(idx, \(.) {
    z <- assay(y)[, .]
    z <- as.matrix(z)
    z[z < 0] <- 0
    prop.table(z, 1) # equivalent to x/sum(x)
  })
  # compute entropy...
  h <- \(p) {
    p <- p[p > 0]
    n <- log(length(p), base = base)
    -sum(p*log(p, base = base))/n
  }
  # ...across clusters for every sample
  hs <- sapply(p, apply, 1, h)
  z <- aggregateAcrossCells(x, 
    x[[j]], coldata.merge = FALSE,
    use.assay.type = assay, statistics = fun)
  res <- data.frame(
    row.names = NULL, entropy = c(hs),
    marker_id = rep(rownames(hs), ncol(hs)),
    sample_id = rep(colnames(hs), each = nrow(hs)))
  res[is.na(res)] <- 0
  return(res)
}

.score <- \(se, da, ds) {
  ms <- unique(se$marker_id, ds$marker_id)
  ss <- c(by(ds, ds$marker_id, \(.) mean(abs(.$logFC))))
  ts <- c(by(se, se$marker_id, \(.) mean(1-.$entropy)))
  df <- data.frame(
    row.names = NULL, marker_id = ms,
    state_score = ss[ms], type_score = ts[ms])
}

.ei <- \(x, 
  sample = "sample_id",
  group = "condition")
{
  ids <- unique(x[[sample]])
  idx <- match(ids, x[[sample]])
  ei <- colData(x)[idx, c(sample, group)]
  ei <- data.frame(ei, row.names = NULL)
} 
# setup experiment info 

.da <- \(x, 
  cluster = names(cluster_codes(x))[1],
  sample = "sample_id", group = "condition")
{
  # setup model & contrast matrix
  ei <- .ei(x, sample, group)
  dm <- createDesignMatrix(ei, group)
  cm <- createContrast(c(0, 1))
  # run differential abundance (DA) analysis
  se <- diffcyt(x, ei,
    clustering_to_use = cluster, 
    design = dm, contrast = cm,  verbose = FALSE, 
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR")$res
  data.frame(rowData(se))
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

.calculate_wilcox_exprs_score <- \(x, clustering_to_use = names(cluster_codes(sce)[1]), 
                                  assay_to_use = "exprs", fun = "mean", penalty = T){
  # wilcox sum rank p-value
  clusters <- levels(cluster_ids(x, clustering_to_use))
  exprs_mm <- assay(x, assay_to_use)
  n_genes <- 1:nrow(x)
  wilcox_per_cluster <- sapply(clusters, \(j){
    temp <- x
    idx <- which(colData(temp)$cluster_id == j)
    colData(temp)[idx,]$cluster_id <- 1
    colData(temp)[-idx,]$cluster_id <- 2
    wilcox_res <- sapply(n_genes, function(i){
      wilcox.test(exprs_mm[i,] ~ colData(temp)$cluster_id)$p.value
    })
  })
  rownames(wilcox_per_cluster) <- rownames(x)
  mean_wilcox <- apply(wilcox_per_cluster, 1, mean, na.rm=TRUE)
  
  
  # expression value
  if(penalty == T){
    colData(x)[,clustering_to_use] <- cluster_ids(x, clustering_to_use)
    pb <- muscat::aggregateData(x, by = clustering_to_use, assay = assay_to_use, fun = "mean")
    minmax_norm <- \(x) {
      (x - min(x)) / (max(x) - min(x))
    }
    norm_mat <- apply(assay(pb), 2, minmax_norm)
    max_mat <- apply(norm_mat, 1, max)
    score <- max_mat*(-log(mean_wilcox))
    
  }else{
    score <- -log(mean_wilcox)
    
  }
  return(score)
}


.one_vs_rest_logFC <- \(x, clustering_to_use = names(cluster_codes(x)[1]), 
                        assay_to_use = "exprs", fun = "mean"){
  colData(x)[,clustering_to_use] <- cluster_ids(x, clustering_to_use)
  pb <- aggregateData(x, by = clustering_to_use, assay = assay_to_use, fun = fun)
  clusters <- levels(cluster_ids(x, clustering_to_use))
  exprs_mm <- assay(x, assay_to_use)
  n_genes <- 1:nrow(x)
  custom_function <- \(x) {
    \(i) x[i] / mean(x[x != x[i]])
  }
  logFC_per_cluster <- t(apply(assay(pb), 1, function(x) sapply(seq_along(x), custom_function(x))))
  return(logFC_per_cluster)
}

.Silhouette_score <- \(x, dr = "PCA", feature = "type", 
                       clustering_to_use = "meta20", 
                       dim = 10, seed=1234,
                       fun = mean){
  x <- runDR(x, dr = dr, feature = feature)
  x <- cluster(x, 
                 xdim = dim, ydim = dim,
                 features = feature, 
                 seed = seed, verbose = FALSE)
  colData(x)[,clustering_to_use] <- cluster_ids(x, clustering_to_use)
  sil.approx <- approxSilhouette(reducedDim(x, dr), clusters=colData(x)[,clustering_to_use])
  sil.data <- as.data.frame(sil.approx)
  sil.data$closest <- factor(ifelse(sil.data$width > 0, sil.data$cluster, sil.data$other))
  fun(sil.data$width)
}


