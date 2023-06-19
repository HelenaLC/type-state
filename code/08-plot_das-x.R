args <- list(
    list.files("outs", "^das.*", full.names = TRUE),
    "plts/das-x.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(.) .[c("t", "s", "das", "p_val", "p_adj")])

typ <- vapply(res, \(.) .$das[1], character(1))
ds <- !(da <- grepl("^DA", typ))
da <- do.call(rbind, res[da])
ds <- do.call(rbind, res[ds])

ggplot(ds, aes(p_val, col = das)) +
    facet_grid(s ~ t) +
    geom_density()
