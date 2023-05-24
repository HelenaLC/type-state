suppressPackageStartupMessages(
  {
   library(bluster)
  }
)

fun <- \(x, 
                    sample_id = "sample_id", 
                    cluster_id = "cluster_id", 
                    fun = mean, dr = "PCA"){
  sil <- sapply(unique(colData(x)[,sample_id]), \(i){
    idx <- which(colData(x)[,sample_id] == i)
    temp <- x[, idx]
    sil.approx <- approxSilhouette(reducedDim(temp, dr), clusters=colData(temp)[,cluster_id])
    sil.data <- as.data.frame(sil.approx)
    sil.data$closest <- factor(ifelse(sil.data$width > 0, sil.data$cluster, sil.data$other))
    fun(sil.data$width)
  })
  return(fun(sil))
}