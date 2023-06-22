# args <- list(
#     list.files("outs", "^das.*", full.names = TRUE),
#     "plts/das-x.pdf")

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

plt <- ggplot(ds, aes(p_val, col = das)) +
    facet_grid(s ~ t, scales = "free", labeller = \(.) label_both(.)) +
    geom_density() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(args[[2]], plt, units = "cm", width = 20, height = 15)
