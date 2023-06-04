####### load required packages #######
suppressPackageStartupMessages({
  library("splatter")
})

fun <- \()
{
  set.seed(seed <- 42)
  vcf <- mockVCF()
  gff <- mockGFF(n.genes = 2000, seed = seed)
  params.group <- newSplatPopParams(batchCells = 500,
                                    similarity.scale = 0,
                                    de.prob = 0.3,
                                    de.facLoc = 1, 
                                    de.facScale = 1,
                                    de.downProb = 0,
                                    group.prob = c(1/3, 1/3, 1/3),
                                    condition.prob = c(0.5, 0.5),
                                    cde.prob = 0,
                                    cde.facLoc = 0, 
                                    cde.facScale = 0)
  
  x <- splatPopSimulate(vcf = vcf, gff = gff, 
                        params = params.group, 
                        sparsify = FALSE)
 
  return(x)
}

