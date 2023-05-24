suppressPackageStartupMessages({
  library("splatter")
})

fun <- \()
{
  params.group <- newSplatPopParams(batchCells = 50,
                                    similarity.scale = 0,
                                    de.prob = 0.5,
                                    de.facLoc = 0.5, 
                                    de.facScale = 0.5,
                                    group.prob = c(1/3, 1/3, 1/3),
                                    condition.prob = c(0.5, 0.5),
                                    cde.prob = 0.5, 
                                    cde.facLoc = 0.5, 
                                    cde.facScale = 0.5)
  
  sim.state_type <- splatPopSimulate(vcf = vcf, gff = gff, params = params.group, 
                                     sparsify = FALSE)
  
  sim.state_type <- logNormCounts(sim.state_type)
  sim.state_type <- runPCA(sim.state_type, ncomponents = 10)
  colData(sim.state) <- DataFrame(sample_id = sim.state_type$Sample,
                                  cluster = sim.state_type$Group,
                                  condition = sim.state_type$Condition,
                                  row.names = NULL)
}