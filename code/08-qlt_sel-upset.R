# args <- list(
#     list.files("outs", "^sel-", full.names=TRUE),
#     "plts/sel-upset.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ComplexUpset)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)
p <- lapply(split(df, df$dat), \(fd) {
    gs <- lapply(
        split(df, df$sel), 
        \(.) .[.$sel_val,"gene_id"])
    tf <- UpSetR::fromList(gs)
    ComplexUpset::upset(tf,
        intersect = colnames(tf),
        set_sizes=FALSE, n_intersections=20, 
        height_ratio=1, min_degree=1, min_size=10, 
        name=df$dat[1]) 
    
}) |> wrap_plots()


ggsave(args[[2]], p, width=15, height=10, units="cm")