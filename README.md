### setup

- workflow was implemented and last executed successfully with<br>
  **R v4.4.1 with Bioc 3.20, and Python v3.11.3 with Snakemake v7.26.0**
- R version and library have to be specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `.Rprofile` is used for handling and printing command line arguments
- `logs/` capture `.Rout` files from `R CMD BATCH` executions
- `data/` contains any synthetic and real data
- intermediate results are generated in `outs/` 
- visualizations are generated in `plts/`

### workflow

- `<x>` denotes a wildcard, namely: `t`ype, `s`tate, `b`atch,  
  `sim`ulation, `sco`re, `sel`ection, `sta`tistic,   
  `das` = differential state analysis method

- `00-get_sim/dat.R`
  - **out:** for simulations, `data/sim/00-raw/t<t>,s<s>,b<b>.rds`,<br>
    for real data, `data/dat/00-raw/<did>.rds` (`<did>` = dataset identifier)
  - synthetic data generation (`splatter::splatPopSimulate()`)
  - hereafter, `t<t>,s<s>,b<b>` = `<sim>`

- `01-pro_sim/dat.R`
  - **in:** `data/sim|dat/00-raw/<sim|dat>.rds`<br>
    **out:** `data/sim|dat/01-fil/<sim|dat>.rds`
  - minimal filtering keeping genes with count > 1  
    in ≥ 10 cells, and cells with ≥ 10 detected genes
  - log-library size normalization (`scater::logNormCounts()`)
  - highly variable gene (HVG) selection (`scran::modelGeneVar()`)
  - principal component analysis (PCA) using HVGs (`scater::runPCA()`)

- `02-sco.R`
  - **in:** `data/sim|dat/01-fil/<sim|dat>.rds`<br>
    **out:** `outs/sim|dat/sco-<sim|dat>,<sco>.rds`
  - source method from one of `02-sco-<sco>.R`
  - compute gene-level metrics to quantify type-/state-specificity 

- `03-sel.R`
  - **in:** `outs/sim|dat/sco-<sim|dat>,<sco>.rds`<br>
    **out:** `outs/sim|dat/sel-<sim|dat>,<sco>.rds`
  - source method from one of `03-sel-<sel>.R`
  - select genes for reprocessing

- `04-rep.R`
  - **in:** `outs/sim|dat/sco-<sim|dat>,<sco>.rds`<br>
    **out:** `data/sim|dat/02-rep/<sim|dat>,<sel>.rds`
  - data reprocessing (PCA, clustering, reduction)
  
- `05-sta.R`
  - **in:** `data/sim|dat/02-rep/<sim|dat>,<sel>.rds`<br>
    **out:** `outs/sim|dat/sta-<sim|dat>,<sel>,<sta>.rds`
  - source method from on of `05-sta-<sta>.R`
  - compute evaluation statistics

- `06-das.R`
  - **in:** `data/sim|dat/02-rep/<sim|dat>,<sel>.rds`<br>
    **out:** `outs/sim|dat/das-<sim|dat>,<sel>,<das>.rds`
  - source method from one of `06-das-<das>.R`
  - perform differential state analysis (DSA)

- `07-eva.R`
  - standalone script applied to experimental data only
  - collects results across all feature selection strategies,<br>
    selects [10, 20, ..., 90\%] for top-rank features, and recomputes<br>
    evaluation statistics for accordingly reprocessed data (PCA, clustering)

- `08-plt_<out>-<plt>.R`
  - **in:** `outs/sim/<out>.rds`<br>
    **out:** `plt/sim/<out>-<plt>.pdf`
  - e.g., `08-plt_das-F1.pdf` collects all DSA results<br>
    (`outs/sim/das-<sim>,<sel>,<das>.rds`) and plots F1 scores
  - visualization of synthetic data analysis results

- `08-qlt_<out>-<qlt>.R`
  - **in:** `outs/dat/<out>.rds`<br>
    **out:** `plt/dat/<out>-<qlt>.pdf`
  - visualization of experimental data analysis results
  
- `09-aes.R`
  - sourced to fix the order of feature scores (`SCO`),<br>
    ground truth-based (`DES`) and other selections (`SEL`),<br>
    and differential state analysis methods (`DAS`) across plots

- `10-session_info.R`
  - lists and may be used to install all R packages used<br>
    (across CRAN, GitHub, and Bioconductor), and writes the<br>
    corresponding `sessionInfo()` output to `session_info.txt`