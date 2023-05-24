suppressPackageStartupMessages(
  {
    library(limma)
    library(muscat)
    library(SingleCellExperiment)
  }
)

fun <- \(x, 
         sample = "sample_id", 
         cluster = "cluster_id",
         assay = "counts",
         condition = "condition",
         fun = "sum"){
  
  res <- lapply(unique(colData(x)[,cluster]), \(i){
    temp <- x[, which(colData(x)[,cluster] == i)]
    sce <- SingleCellExperiment(assays=list(counts=assay(temp,
                                                         assay)),
                                colData=DataFrame(sample_id = colData(temp)[,sample]))
    
    pb <- aggregateData(sce, assay = assay, fun=fun, by=c("sample_id"))
    ei <- .ei(temp, sample, condition)
    dm <- createDesignMatrix(ei, condition)
    cm <- createContrast(c(0, 1))
    dge <- DGEList(counts = assay(pb), remove.zeros = T)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design = dm, plot = F)
    fit <- lmFit(v, design = dm)
    fit2 <- contrasts.fit(fit, contrasts = cm)
    fit2 <- eBayes(fit2)
    limma_out <- topTable(fit2, number = nrow(dge))
  })
  names(res) <- unique(colData(x)[,cluster])
  final <- sapply(res, \(ds){
    idx <- match(rownames(x), rownames(ds))
    ds$logFC[idx]
  })
  rownames(final) <- rownames(x)
  rowMeans(abs(final))
}