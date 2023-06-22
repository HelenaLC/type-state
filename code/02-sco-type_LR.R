suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(lmtest)
})

# Code adapted from Seurat
fun <- \(x, 
    aggregate_to_use = mean, 
    assay_to_use = "logcounts", 
    cluster_to_use = "cluster_id") {

    y <- assay(x, "logcounts")
    ids <- unique(x$cluster_hi)
    res <- sapply(ids, \(k) {
        tmp <- x
        id <- tmp$cluster_hi
        j <- !(i <- id == k)
        ij <- c(which(i), which(j))
        df <- data.frame(row.names = ij)
        df[i, "group"] <- "Group1"
        df[j, "group"] <- "Group2"
        df$group <- factor(df$group)
        z <- y[, as.numeric(rownames(df)), drop = FALSE]
        p_val <- sapply(seq_len(nrow(x)), \(gene){
            model.data <- cbind(GENE = z[gene, ], df)
            fmla <- as.formula(object = "group ~ GENE")
            fmla2 <- as.formula(object = "group ~ 1")
            model1 <- glm(formula = fmla, 
                data = model.data, 
                family = "binomial")
            
            model2 <- glm(formula = fmla2, 
                data = model.data, 
                family = "binomial")
            
            lrtest <- lrtest(model1, model2)
            lrtest$Pr[2]
        })
    })
    #rownames(res) <- rownames(x)
    res <- -log(res)
    apply(res, 1, FUN = mean, na.rm = TRUE)
    
}