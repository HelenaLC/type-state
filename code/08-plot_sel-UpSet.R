suppressPackageStartupMessages({
    library(UpSetR)
    library(ggplot2)
    library(dplyr)
    library(ggplotify)
    library(cowplot)
})


sel <- lapply(args[[1]], readRDS)
names(sel) <- sapply(sel, \(.) .$sel[1])
df <- do.call(rbind, sel) 


plt <- list()
for (t in unique(df$t)) {
    for(s in unique(df$s)){
        tmp <- df[df$t == t & df$s == s,]
        lst <- split(tmp, tmp$sel)
        sel_idx <- lapply(lst, \(x) x$gene_id[x$sel_val == TRUE])
        plt[[paste0("t",t,",s",s)]] <-
            as.ggplot(upset(fromList(sel_idx),
                order.by = "freq",
                nsets = length(sel_idx),
                set_size.show = TRUE,
                scale.sets = "identity"))
        
    }
}


res <- plot_grid(plotlist = plt, labels = names(plt), ncol = 6, nrow = 6)
ggsave(filename = args[[2]], plot = res, width = 40, height = 30)








