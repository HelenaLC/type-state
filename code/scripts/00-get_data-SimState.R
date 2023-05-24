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
                                    de.prob = 0,
                                    de.facLoc = 0, 
                                    de.facScale = 0,
                                    group.prob = c(1/3, 1/3, 1/3),
                                    condition.prob = c(0.5, 0.5),
                                    cde.prob= 0.9,
                                    cde.facLoc = 1, 
                                    cde.facScale = 1)
  
  sim.state <- splatPopSimulate(vcf = vcf, gff = gff, params = params.group, 
                                sparsify = FALSE)
  
  # sim.state <- logNormCounts(sim.state)
  # sim.state <- runPCA(sim.state, ncomponents = 10)
  
  colData(sim.state) <- DataFrame(sample_id = sim.state$Sample,
                                  cluster = sim.state$Group,
                                  condition = sim.state$Condition,
                                  row.names = NULL)
  return(sim.state)                                  
}
