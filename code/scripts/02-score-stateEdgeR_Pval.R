suppressPackageStartupMessages(
    {
        library(limma)
        library(muscat)
        library(SingleCellExperiment)
        source("code/scripts/utils.R")
        library(edgeR)
    }
)

fun <- \(x, 
         sample = "sample_id", 
         cluster = "cluster_id",
         assay_to_use = "counts",
         condition = "condition",
         fun = "sum",
         min_cells = 3,
         min_samples = 1){
    
    counts <- assay(x, "counts")
    tf <- counts >= min_cells
    ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
    x <- x[ix_keep, ]
    
    res <- lapply(unique(colData(x)[,cluster]), \(i) {
        temp <- x[, which(colData(x)[,cluster] == i)]
        sce <- SingleCellExperiment(assays = list(counts = assay(temp, assay_to_use)),
                                    colData = DataFrame(sample_id = colData(temp)[, sample],
                                                        condition = colData(temp)[, condition]))
        # drop samples without any cells
        sce$sample_id <- droplevels(sce$sample_id)
        # split cell indices by sample
        idx <- split(seq(ncol(sce)), sce$sample_id)
        
        # construct design matrix
        sids <- unique(sce$sample_id)
        idx <- match(sids, sce$sample_id)
        ei <- data.frame(
            row.names = NULL, 
            sample_id = sids,
            condition = sce$condition[idx])
        group <- factor(ei$condition)
        # if more than two levels in design matrix, do DS analysis
        if (length(levels(group)) > 1) {
            colData(sce)$pb_group <- paste0(colData(sce)$sample_id, 
                                            "-",
                                            colData(sce)$condition)
            sce_counts <- assay(sce, "counts")
            pb_counts <- t(rowsum(t(sce_counts), colData(sce)$pb_group))
            pb_samples <- colnames(pb_counts)
            pb_split <- do.call(rbind, strsplit(pb_samples, "-"))
            group <- pb_split[, 2]
            samples <- pb_split[, 1]
            y <- DGEList(counts = pb_counts, 
                         group = group,
                         samples = samples)
            #keep <- filterByExpr(y)
            #y <- y[keep, , keep.lib.sizes=TRUE]
            y <- calcNormFactors(y)
            design <- model.matrix(~group)
            y <- estimateDisp(y, design)
            fit <- glmFit(y, design)
            lrt <- glmLRT(fit)
            
            return(lrt$table)
        } else {
            return(NULL)
        }
        
    })
    
    lst <- lapply(res, \(ds) {
        if (!is.null(ds)) {
            idx <- match(rownames(ds), rownames(x))
            ss <- data.frame(row.names = rownames(x),
                             score = replicate(nrow(x), 0))
            ss$score[idx] <- 1 - ds$PValue
            return(ss)
            
        } 
    }) 
    final <- do.call(cbind, lst)

    return(rowMeans(abs(final)))
    
    
}
