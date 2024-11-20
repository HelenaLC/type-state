### setup

- R version and library have to be specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `logs` captures `.Rout` files from `R CMD BATCH` executions
- `data` contains any synthetic and real data
- intermediate results are generated in `outs` 
- visualizations are generated in `plts`

### workflow

- `<x>` denotes a wildcard, namely: `t`ype, `s`tate, `b`atch,  
  `sim`ulation, `sco`re, `sel`ection, `sta`tistic,   
  `das` = differential state analysis method

- `00-sim_data.R`
  - **out:** `data/00-sim/t<t>,s<s>,b<b>.rds`<br>
    (hereafter, `t<t>,s<s>,b<b>` = `<sim>`)
  - synthetic data generation (`splatter::splatPopSimulate()`)

- `01-fil_data.R`
  - **in:** `data/00-sim/<sim>.rds`<br>
    **out:** `data/01-fil/<sim>.rds`
  - minimal filtering keeping genes with count > 1  
    in ≥ 10 cells, and cells with ≥ 10 detected genes
  - log-library size normalization (`scater::logNormCounts()`)
  - highly variable gene (HVG) selection (`scran::modelGeneVar()`)
  - principal component analysis (PCA) using HVGs (`scater::runPCA()`)

- `02-sco.R`
  - **in:** `data/01-fil/<sim>.rds`<br>
    **out:** `outs/sco-sim-<sim>,<sco>.rds`
  - source method from one of `02-sco-<sco>.R`
  - compute gene-level metrics to quantify type-/state-specificity 

- `03-sel.R`
  - **in:** `outs/sco-sim-<sim>,<sco>.rds`<br>
    **out:** `outs/sel-sim-<sim>,<sco>.rds`
  - source method from one of `03-sel-<sel>.R`
  - select genes for reprocessing

- `04-rep.R`
  - **in:** `outs/sco-sim-<sim>,<sco>.rds`<br>
    **out:** `data/02-rep/sim-<sim>,<sel>.rds`
  - data reprocessing (PCA, clustering, reduction)
  
- `05-sta.R`
  - **in:** `data/02-rep/sim-<sim>,<sel>.rds`<br>
    **out:** `outs/sta-sim-<sim>,<sel>,<sta>.rds`
  - source method from on of `05-sta-<sta>.R`
  - compute evaluation statistics

- `06-das.R`
  - **in:** `data/02-rep/sim-<sim>,<sel>.rds`<br>
    **out:** `outs/das-sim-<sim>,<sel>,<das>.rds`
  - source method from one of `06-das-<das>.R`
  - perform differential abundance/state analysis 
