# retrieve dataset
source(args[[1]])
x <- fun()

# number of features & observations
dim(x)

s <- !is.null(x$sample_id)
k <- !is.null(x$cluster_id)
if (b) {
  table(x$sample_id)
  if (k) {
    table(x$cluster_id)
    table(x$sample_id, x$cluster_id)
  }
} else if (k) {
  table(x$cluster)
}

dimnames(x) <- list(
  paste0("gene", seq_len(nrow(x))),
  paste0("cell", seq_len(ncol(x))))

saveRDS(x, args[[2]])