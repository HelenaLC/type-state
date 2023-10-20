suppressPackageStartupMessages({
    library(splatter)
})

# setup
set.seed(seed <- 7043)
t <- as.numeric(wcs$t)/100
s <- as.numeric(wcs$s)/100
b <- as.numeric(wcs$b)/100

p <- newSplatPopParams(
    de.facLoc=t,
    de.facScale=0.001,
    cde.facLoc=s,
    cde.facScale=0.001,
    batch.facLoc=b,
    bcv.common=1.5,
    similarity.scale=10,
    de.prob=0.2, 
    cde.prob=0.2, 
    batchCells=400,
    group.prob=rep(1/3, 3),
    condition.prob=c(0.5, 0.5))

vcf <- mockVCF(n.samples=6, seed=seed)
gff <- mockGFF(n.genes=2e3, seed=seed)

# data generation
x <- splatPopSimulate(
    params=p, vcf=vcf, gff=gff,
    verbose=FALSE, sparsify=FALSE)

# standardize cell metadata
colData(x) <- DataFrame(
    cluster_id=x$Group,
    sample_id=x$Sample,
    group_id=x$Condition)
for (. in names(colData(x)))
    x[[.]] <- factor(x[[.]])

# store gene/cell identifiers
.f <- \(x, .) {
    switch(., row={lab="gene"; dim=nrow(x)}, col={lab="cell"; dim=ncol(x)})
    sprintf("%s%s", lab, formatC(x=seq_len(dim), width=nchar(dim), flag=0))
}
rowData(x)$gene_id <- rownames(x) <- .f(x, "row")
colData(x)$cell_id <- colnames(x) <- .f(x, "col")

# store simulation parameters
metadata(x) <- list(t=t, s=s, b=b)

# save data
saveRDS(x, args$res)