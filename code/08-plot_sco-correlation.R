#args <- list(list.files("outs", "^sco-.*", full.names = TRUE), "plts/sco-density.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
    library(stringr)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) 

## define type marker genes
gde <- df[grep("GroupDE", names(df))]
mg <- sapply(seq_len(ncol(gde)), \(i){
    not_i <- setdiff(seq_len(ncol(gde)), i)
    mgk <- sapply(not_i, \(j) {
        log(gde[,i]/gde[,j], base = 2)
    })
    rowMeans(mgk)
})

df$mg <- apply(abs(mg), 1, mean)

## define state genes
cde <- df[grep("ConditionDE", names(df))]
cg <- sapply(seq_len(ncol(cde)), \(i){
    not_i <- setdiff(seq_len(ncol(cde)), i)
    cgk <- sapply(not_i, \(j) {
        log(cde[,i]/cde[,j], base = 2)
    })
    rowMeans(cgk)
})
df$cg <- apply(abs(cg), 1, mean)

state <- df[str_detect(df$sco, "state"),]
lst <- split(df, list(df$sco, df$t, df$s, df$b))
scr <- lapply(lst, \(x) {
    if (str_detect(x$sco[1], "state")) {
        cor <- cor(x$sco_val, x$cg, method = "spearman")
    } else {
        cor <- cor(x$sco_val, x$mg, method = "spearman")
    }
    data.frame(sco = x$sco[1], t = x$t[1], 
        s = x$s[1], b = x$b[1], cor = cor)
}) 
res <- do.call(rbind, scr)

gg <- ggplot(res) +
    facet_grid( ~ sco) +
    geom_tile(col = "white", aes(t, s, fill = cor)) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    labs(x = "type effect", y = "state effect") +
    scale_fill_distiller(NULL,
        palette = "RdBu", na.value = "lightgrey",
        limits = c(-1, 1), n.breaks = 3, direction = 1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA))

ggsave(args[[2]], gg, units = "cm", width = 30, height = 25)


