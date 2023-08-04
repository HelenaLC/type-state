.Fstat <- \(x) {
    y <- assay(x, "logcounts")
    cd <- data.frame(colData(x))
    f <- ~ sample_id + cluster_hi
    mm <- model.matrix(f, data = cd)
    rownames(mm) <- colnames(x)
    res <- apply(y, 1, \(g) {
        fit <- lmFit(g, mm)
        fit <- eBayes(fit, trend = FALSE)
        cs <- colnames(fit$cov.coefficients)
        nan <- !colnames(mm) %in% cs
        if (any(nan)) {
            fit <- lmFit(g, mm[, !nan])
            fit <- eBayes(fit, trend = FALSE)
        }
        cs <- grep("cluster", cs)
        topTable(fit, coef = cs, sort.by = "none")$F
    })
    return(res)
}

.edgeR <- \(x) {
    ids <- colData(x)[c("sample_id", "cluster_lo")]
    y <- aggregateAcrossCells(x, ids)
    idx <- split(seq(ncol(y)), y$cluster_lo)
    res <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$group_id)
        if (length(gs) == 2) {
            mm <- model.matrix(~z$group_id)
            z <- DGEList(assay(z))
            z <- calcNormFactors(z)
            z <- estimateDisp(z, mm)
            fit <- glmQLFit(z, mm)
            lrt <- glmQLFTest(fit, 2)
            tbl <- topTags(lrt, n = Inf, sort.by = "none")$table
            tbl <- rename(tbl, p_val = "PValue", p_adj = "FDR")
            data.frame(
                row.names = NULL,
                gene = rownames(tbl), 
                cluster_id = k, tbl)
        }
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) {
        # average across clusters
        res <- group_by(res, gene)
        res <- summarize(res, mean(-log(p_adj)))
        out <- numeric(nrow(x))
        names(out) <- rownames(x)
        out[res$gene] <- res[[2]]
        return(out)
    }
}

.limma <- \(x) {
    # pseudo-bulks by cluster-sample
    y <- aggregateAcrossCells(x,
        ids = colData(x)[c("cluster_lo", "sample_id")],
        use.assay.type = "logcounts", statistics = "mean")
    # test for DS by cluster
    idx <- split(seq(ncol(y)), y$cluster_lo)
    tbl <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        z$sample_id <- droplevels(z$sample_id)
        ns <- table(z$sample_id, z$group_id)
        if (!any(colSums(ns) == 0)) {
            mm <- model.matrix(~ z$group_id)
            ids <- levels(z$sample_id)
            rownames(mm) <- ids
            ws <- table(z$sample_id)[ids]
            fit <- lmFit(assay(z), design = mm, weights = ws)
            fit <- contrasts.fit(fit, c(0, 1))
            fit <- eBayes(fit, trend = FALSE)
            tbl <- topTable(fit, sort = "none", n = Inf)
            rnm <- match(c("P.Value", "adj.P.Val"), names(tbl))
            names(tbl)[rnm] <- c("p_val", "p_adj")
            data.frame(
                cluster_id = k,
                gene = rownames(tbl), 
                tbl, row.names = NULL)
        }
    })
    tbl <- do.call(rbind, tbl)
    if (is.null(tbl)) return(NULL)
    # average across clusters
    avg <- tbl %>% 
        group_by(gene) %>%
        summarize(mean(-log(p_adj)))
    res <- numeric(nrow(x))
    names(res) <- rownames(x)
    res[avg$gene] <- avg[[2]]
    return(res)
}


.entropy <- \(x) {
    # compute pseudo-bulks by cluster-sample
    x <- SingleCellExperiment(assays = list(counts = assay(x, "counts")),
        colData = DataFrame(sample_id = x$sample_id,
            condition = x$group_id,
            cluster_id = x$cluster_hi))
    
    y <- aggregateAcrossCells(x, 
        colData(x)[c("cluster_id", "sample_id")], 
        coldata.merge = FALSE,
        use.assay.type = "counts", 
        statistics = "sum")
    
    # split cell indices by samples
    idx <- split(seq(ncol(y)), y$sample_id)
    # compute proportions across clusters
    # (matrix of dim. features x clusters)
    p <- lapply(idx, \(.) {
        z <- assay(y)[, .]
        z <- as.matrix(z)
        z[z < 0] <- 0
        prop.table(z, 1) # equivalent to x/sum(x)
    })
    # compute entropy...
    h <- \(p) {
        p <- p[p > 0]
        n <- log(length(p), base = 2)
        -sum(p*log(p, base = 2))/n
    }
    # ...across clusters for every sample
    hs <- sapply(p, apply, 1, h)
    z <- aggregateAcrossCells(x, 
        x$sample_id, coldata.merge = FALSE,
        use.assay.type = "counts", statistics = "sum")
    res <- data.frame(
        row.names = NULL, entropy = c(hs),
        marker_id = rep(rownames(hs), ncol(hs)),
        sample_id = rep(colnames(hs), each = nrow(hs)))
    # when cluster size is small,
    # a gene is likely to only express in one cluster
    res[is.na(res)] <- 0
    c(by(res, res$marker_id, \(.) mean(1-.$entropy)))
}

.DA_edgeR <- \(x) {
    y <- table(x$sample_id, x$cluster_re)
    y <- t(as.matrix(unclass(y)))
    
    ids <- colnames(y)
    idx <- match(ids, x$sample_id)
    df <- data.frame(
        row.names = ids,
        sample_id = ids,
        group_id = x$group_id[idx])
    mm <- model.matrix(~ group_id, df)
    dge <- DGEList(y, samples = df)
    dge <- estimateDisp(dge, mm, trend = "none")
    tbl <- tryCatch({
        fit <- glmQLFit(dge, mm, robust = TRUE, abundance.trend = FALSE)
        glmQLFTest(fit, coef = ncol(mm))},
        error = function(e) NULL)
    res <- if (!is.null(tbl)) {
        tbl <- data.frame(topTags(tbl, sort.by = "none"))
        idx <- match(c("PValue", "FDR"), colnames(tbl))
        colnames(tbl)[idx] <- c("p_val", "p_adj")
        cell_n <- tabulate(x$cluster_re)
        cell_i <- split(seq(ncol(x)), x$cluster_re)
        DataFrame(
            row.names = NULL, tbl,
            cluster_re = rownames(tbl), 
            cell_n[rownames(tbl)], 
            cell_i = I(cell_i)[rownames(tbl)])
    }
    res <- data.frame(res, DAS = "DA_edgeR")
    return(res)
}


.DA_limma <- \(x) {
    y <- table(x$sample_id, x$cluster_re)
    y <- t(as.matrix(unclass(y)))
    
    ids <- colnames(y)
    idx <- match(ids, x$sample_id)
    df <- data.frame(
        row.names = ids,
        sample_id = ids,
        group_id = x$group_id[idx])
    mm <- model.matrix(~ group_id, df)
    
    y <- tryCatch(
        voom(DGEList(y), mm, plot = FALSE),
        error = function(e) NULL)
    df <- if (!is.null(y)) {
        fit <- lmFit(y, mm)
        fit <- contrasts.fit(fit, c(0, 1))
        fit <- eBayes(fit, trend = FALSE)
        
        tbl <- topTable(fit, number = Inf, sort.by = "none")
        idx <- match(c("P.Value", "adj.P.Val"), names(tbl))
        names(tbl)[idx] <- c("p_val", "p_adj")
        
        cell_n <- tabulate(x$cluster_re)
        cell_i <- split(seq(ncol(x)), x$cluster_re)
        DataFrame(
            row.names = NULL, tbl,
            cluster_re = rownames(y), 
            cell_n, cell_i = I(cell_i))
    }
    df <- data.frame(df, DAS = "DA_limma")
    return(df)
}

.DA_milo <- \(x) {
    cd <- data.frame(colData(x))
    df <- distinct(cd[c("sample_id", "group_id")])
    rownames(df) <- df$sample_id
    
    y <- makeNhoods(buildGraph(Milo(x), d = ncol(reducedDim(x, "PCA"))))
    y <- countCells(y, "sample_id", cd)
    y <- calcNhoodDistance(y, d = ncol(reducedDim(x, "PCA")))
    res <- testNhoods(y, ~ group_id, df)
    
    idx <- match(c("PValue", "FDR"), names(res))
    names(res)[idx] <- c("p_val", "p_adj")
    res <- data.frame(res, DAS = "DA_milo")
    return(res)
}


.DS_edgeR <- \(x) {
    ids <- colData(x)[c("sample_id", "cluster_re")]
    y <- aggregateAcrossCells(x, ids)
    idx <- split(seq(ncol(y)), y$cluster_re)
    res <- lapply(names(idx), \(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$group_id)
        if (length(gs) == 2 & length(idx[[k]]) > 2) {
            mm <- model.matrix(~ z$group_id)
            z <- DGEList(assay(z))
            z <- calcNormFactors(z)
            z <- estimateDisp(z, mm)
            fit <- glmQLFit(z, mm)
            lrt <- glmQLFTest(fit)
            tbl <- topTags(lrt, n = Inf, sort.by = "none")
            data.frame(
                row.names = NULL,
                gene = rownames(tbl), 
                cluster_re = k, tbl)
        }
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) {
        idx <- match(c("PValue", "FDR"), names(res))
        names(res)[idx] <- c("p_val", "p_adj")
    }
    res <- data.frame(res, DAS = "DS_edgeR")
    return(res)
}

.DS_lemur <- \(x) {
    fit <- lemur(x, design = ~ group_id)
    fit <- align_harmony(fit) 
    cd1 <- unique(x$group_id)[1]
    cd2 <- unique(x$group_id)[2]
    fit <- test_de(fit, 
        contrast = cond(group_id = cd1) - cond(group_id = cd2))
    res <- find_de_neighborhoods(fit, group_by = vars(sample_id, group_id))
    res <- res[res$selection == TRUE,]
    idx <- match(c("pval", "adj_pval"), names(res))
    names(res)[idx] <- c("p_val", "p_adj")
    res <- data.frame(res, DAS = "DS_lemur")
    return(res)
}


.DS_miloDE <- \(x){
    m <- assign_neighbourhoods(x, reducedDim_name = "PCA")
    res <- de_test_neighbourhoods(m, sample_id = "sample_id", 
        design = ~ group_id, covariates = c("group_id"))
    res <- na.omit(res)
    if (!is.null(res)) {
        idx <- match(c("pval", "pval_corrected_across_genes"), names(res))
        names(res)[idx] <- c("p_val", "p_adj")
    }
    res <- data.frame(res, DAS = "DS_miloDE")
    return(res)
}



.sil_k <- \(x) {
    y <- reducedDim(x, "PCA")
    # clusters should be well separated
    # within samples (high value)
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_s <- vapply(idx, \(.) {
        ids <- as.integer(factor(x$cluster_re[.]))
        if (length(unique(ids)) == 1) return(NA)
        res <- silhouette(ids, dist(y[., ]))
        mean(res[, "sil_width"])
    }, numeric(1))
    # samples should be well mixed
    # within clusters (low value)
    res_k <- mean(res_by_s, na.rm = TRUE)

}

.pur_k <- \(x) {
    y <- reducedDim(x, "PCA")
    # cluster purity should be high
    # within samples (high value)
    idx <- split(seq(ncol(x)), x$sample_id)
    res_by_k <- vapply(idx, \(.) {
        ids <- x$cluster_re[.]
        if (length(unique(ids)) == 1) return(NA)
        res <- neighborPurity(y[., ], ids)
        mean(res$purity)
    }, numeric(1))
    
    res_k <- mean(res_by_k, na.rm = TRUE)
}

.FEAST <- \(x) {
    y <- assay(x, "counts")
    y <- as.matrix(y)
    con <- Consensus(Y, k = length(unique(x$cluster_id)))
    res <- cal_F2(y, con$cluster)
    res$F_scores
}

.DUBStepR <- \(x) {
    y <- assay(x, "logcounts")
    colnames(y) <- seq_len(ncol(y))
    #y <- ScaleData(y)
    dub <- DUBStepR(y, optimise.features = FALSE)
    res <- rep(0, nrow(x))
    idx <- match(dub$corr.info$feature.genes, rownames(x))
    res[idx] <- dub$corr.info$corr.range
    return(res)
}

