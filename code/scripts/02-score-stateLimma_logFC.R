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
    assay_to_use = "logcounts",
    condition = "condition",
    fun = "sum",
    min_cells = 3,
    min_samples = 1){
    
    # note: counts are only required for filtering
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
        # compute median expression by marker-sample
        med <- sapply(idx, \(.) rowMedians(assay(sce[, .], "counts")))
        rownames(med) <- rownames(sce)
        
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
            design <- model.matrix(~ group)
            colnames(design) <- levels(group)
            w <- table(sce$sample_id)[sids]
            fit <- lmFit(med, design = design, weights = w)
            fit <- contrasts.fit(fit, c(0, 1))
            fit <- eBayes(fit, trend = TRUE)
            tbl <- limma::topTable(fit, sort = "none", n = Inf)
            rownames(tbl) <- rownames(x)
            return(tbl)
        } else {
            return(NULL)
        }
        
    })
    
    
    lst <- lapply(res, \(ds) {
        if (!is.null(ds)) {
            idx <- match(rownames(ds), rownames(x))
            ss <- data.frame(row.names = rownames(x),
                score = replicate(nrow(x), 0))
            ss$score[idx] <- ds$logFC
            
        } 
    }) 
    final <- do.call(cbind, lst)
    
    rownames(final) <- rownames(x)
    return(rowMeans(final))
    
    
}
