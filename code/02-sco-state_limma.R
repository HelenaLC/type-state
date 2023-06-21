suppressPackageStartupMessages(
    {
        library(limma)
        library(muscat)
        library(SingleCellExperiment)
        library(edgeR)
        library(dplyr)
    }
)

fun <- \(x) {
    
    # note: counts are only required for filtering
    counts <- assay(x, "counts")
    tf <- counts >= 3
    ix_keep <- apply(tf, 1, function(r) sum(r) >= 1)
    x <- x[ix_keep, ]
    
    res <- lapply(unique(x$cluster_lo), \(i) {
        temp <- x[, which(x$cluster_lo == i)]
        sce <- SingleCellExperiment(assays = list(counts = assay(temp, "logcounts")),
            colData = DataFrame(sample_id = temp$sample_id,
                condition = temp$group_id))
        
        # drop samples without any cells
        sce$sample_id <- as.factor(sce$sample_id)
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
            ss$score[idx] <- -log(ds$P.Value)
            
        } 
    }) %>% bind_cols()

    #rownames(final) <- rownames(x)
    rowMeans(lst)
    
        
}
