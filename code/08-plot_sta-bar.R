suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"dat")) readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    group_by(sel, sta, dat) %>%
    summarise_at("sta_val", mean, drop = "none")


plt <- lapply(unique(df$dat), \(x) {
    fd <- df[df$dat == x, ]
    ggplot(fd, aes(x = sel, y = sta_val, col = sel)) +
        geom_bar(stat = "identity", fill = "white") + 
        xlab("Selection method") + ylab("Evaluation values") +
        facet_wrap(~ sta, ncol = 3, scales = "free") + 
        ggtitle(x) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
})
names(plt) <- unique(df$dat)
#ggsave(args[[2]], , units = "cm", width = 30, height = 25)
pdf(args[[2]], width = 12, height = 12, onefile = TRUE)
for (i in names(plt)) {
    print(plt[[i]])
}
dev.off()