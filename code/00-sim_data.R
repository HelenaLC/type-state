suppressPackageStartupMessages({
    library(splatter)
})

set.seed(seed <- 1)

# simulation parameters
t <- as.numeric(wcs$t)/100
s <- as.numeric(wcs$s)/100
b <- as.numeric(wcs$b)/100

p <- newSplatPopParams(
    de.facLoc = t,
    de.facScale = 0.001,
    cde.facLoc = s,
    cde.facScale = 0.001,
    batch.facLoc = b,
    bcv.common = 1.5,
    similarity.scale = 10,
    de.prob = 0.2, cde.prob = 0.2, 
    batchCells = c(400, 400),
    group.prob = rep(1/3, 3),
    condition.prob = c(0.5, 0.5))

vcf <- mockVCF(n.samples = 6, seed = seed)
gff <- mockGFF(n.genes = 2e3, seed = seed)

# data generation
x <- splatPopSimulate(
    params = p, vcf = vcf, gff = gff,
    verbose = FALSE, sparsify = FALSE)

# store simulation parameters
metadata(x) <- list(t = t, s = s, b = b)

# save data
saveRDS(x, args$res)