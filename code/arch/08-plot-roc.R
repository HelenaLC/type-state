suppressPackageStartupMessages({
    library(tidytext)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)

# plt <- ggplot(df, aes(FPR, TPR, col = sco)) + geom_point() + 
#     facet_grid(t ~ s, labeller = \(.) label_both(.)) +
#     theme(axis.text.x = element_text(angle = 45))

a <- ggplot(df, aes(x = factor(t), y = accuracy, col = sco)) +
    geom_boxplot()  +
    xlab("Type effect")

r <- ggplot(df, aes(x =  factor(t), y = TPR, col = sco)) +
    geom_boxplot()  +
    xlab("Type effect")

f <- ggplot(df, aes(x = factor(t), y = FPR, col = sco)) +
    geom_boxplot()  +
    xlab("Type effect") 
    


plt <- a / r / f

ggsave(args[[2]], plt, units = "cm", width = 20, height = 20)