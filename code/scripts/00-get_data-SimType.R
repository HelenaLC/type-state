####### load required packages #######
suppressPackageStartupMessages({
  library("splatter")
})

fun <- \()
{
  set.seed(42)
  vcf <- mockVCF()
  gff <- mockGFF()
  params.group <- newSplatPopParams(batchCells = 50,
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
  
  sim.type <- splatPopSimulate(vcf = vcf, gff = gff, params = params.group, 
                               sparsify = FALSE)
  
  sim.type <- logNormCounts(sim.type)
  sim.type <- runPCA(sim.type, ncomponents = 10)
  
  colData(sim.type) <- DataFrame( sample_id = sim.type$Sample,
                                  cluster = sim.type$Group,
                                  condition = sim.type$Condition,
                                  row.names = NULL)
}

