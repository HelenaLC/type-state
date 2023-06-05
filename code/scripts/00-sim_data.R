suppressPackageStartupMessages({
    library(splatter)
})

t <- as.numeric(wcs$t)/100
s <- as.numeric(wcs$s)/100
b <- as.numeric(wcs$b)/100

set.seed(seed <- 1)
vcf <- mockVCF(n.samples = 6, seed = seed)
gff <- mockGFF(n.genes = 2000, seed = seed)
params.group <- newSplatPopParams(batchCells = c(300, 300),
    batch.facScale = b,
    similarity.scale = 10,
    de.prob = 0.5,
    de.facScale = t,
    group.prob = c(1/3, 1/3, 1/3),
    condition.prob = c(0.5, 0.5),
    cde.prob = 0.5,
    cde.facScale = s,
    bcv.common = 1)

x <- splatPopSimulate(vcf = vcf, gff = gff, 
    params = params.group, 
    sparsify = FALSE)

saveRDS(x, args$res)