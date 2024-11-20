suppressPackageStartupMessages({
  library(lemur)
})

fun <- \(x) {
  set.seed(808)
  ids <- levels(x$group_id)
  fit <- lemur(x, design=~group_id, verbose=FALSE)
  fit <- align_harmony(fit, verbose=FALSE)
  fit <- test_de(fit, cond(group_id=ids[1])-cond(group_id=ids[2]))
  res <- find_de_neighborhoods(fit, vars(sample_id, group_id), verbose=FALSE)
  old <- c("name", "neighborhood", "pval", "adj_pval")
  new <- c("gene_id", "cell_id", "p_val", "p_adj")
  names(res)[match(old, names(res))] <- new
  res$i_nhood <- as.integer(factor(unlist(lapply(res$cell_id, \(.) paste(sort(.), collapse=";")))))
  res$p_nhood <- sapply(res$cell_id, \(.) max(prop.table(table(x[, .]$group_id))))
  res$n_nhood <- length(unique(res$i_nhood))
  return(res)
}