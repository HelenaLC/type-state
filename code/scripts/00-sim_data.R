suppressPackageStartupMessages({
    library(splatter)
})

t <- as.numeric(wcs$t)/100
s <- as.numeric(wcs$s)/100
b <- as.numeric(wcs$b)

set.seed(seed <- 1)
vcf <- mockVCF(n.samples = 6, seed = seed)
gff <- mockGFF(n.genes = 2000, seed = seed)
params.group <- newSplatPopParams(batchCells = 500,
                                  similarity.scale = 10,
                                  de.prob = 0.5,
                                  de.facLoc = t, 
                                  de.facScale = t,
                                  group.prob = c(1/3, 1/3, 1/3),
                                  condition.prob = c(0.5, 0.5),
                                  cde.prob = 0.5,
                                  cde.facLoc = s, 
                                  cde.facScale = s,
                                  bcv.common = b)

x <- splatPopSimulate(vcf = vcf, gff = gff, 
                      params = params.group, 
                      sparsify = FALSE)

saveRDS(x, args$res)