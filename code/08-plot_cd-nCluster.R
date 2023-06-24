suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(reshape2)
    library(SingleCellExperiment)
})

df <- do.call(rbind, lapply(args[[1]], readRDS))
df <- df %>% select(cluster_hi, cluster_lo, t, s)
fd <- melt(df, id.vars = c("t", "s"))
names(fd) <- c("t", "s", "resolution", "id")
gg <- fd %>% 
    group_by(resolution, t, s) %>% 
    summarise(max = max(id))  %>% 
    ggplot(aes(x = resolution, y = max)) + 
    geom_col() +
    facet_grid(t ~ s, labeller = \(.) label_both(.)) +
    geom_text(aes(label = max), vjust = 1, size = 2, colour = "white") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))



ggsave(args[[2]], gg, units = "cm", width = 20, height = 20)
