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

<<<<<<< HEAD

.all_score <- \(x, dim){
  x <- cluster(x, 
                 xdim = dim, ydim = dim,
                 features = rownames(x), 
                 seed = seed, verbose = FALSE)
=======
.all_score <- \(x, dim){
  x <- cluster(x, 
               xdim = dim, ydim = dim,
               features = rownames(x), 
               seed = seed, verbose = FALSE)
>>>>>>> bbc714a (Try different typeness scoring)
  se <- .se(x,  "cluster_id", "sample_id", "exprs", "median")
  se$condition <- x$condition[match(se$sample_id, x$sample_id)]
  da <- .da(x)
  ds <- .ds(x)
  res <- .score(se, da, ds)
  res$marker_class <- marker_classes(x)[res$marker_id]
  return(res)
}

<<<<<<< HEAD
=======


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


>>>>>>> bbc714a (Try different typeness scoring)
