# suppressPackageStartupMessages({
#     library(muscat)
#     library(scran)
#     library(scater)
# })

fun <- \(x) {

  res <- lapply(x, \(m){
    y <- order(m$sco_val, decreasing = TRUE)
  })
  idx <- Reduce(intersect, res)

}
