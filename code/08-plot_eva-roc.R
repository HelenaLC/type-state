suppressPackageStartupMessages({
    library(tidytext)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(caret)
})


res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
roc <- do.call(rbind, res)


p <- ggplot(roc, aes(FPR, TPR, col = sco, shape = sco)) + 
    geom_point() +
    facet_grid(t ~ s, labeller = \(.) label_both(.)) +
    geom_line(aes(linetype = sco))


ggsave(args[[2]], p, units = "cm", width = 40, height = 40)