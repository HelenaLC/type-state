suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
})

res <- lapply(args[[1]], \(x) readRDS(x))
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res) 

p1 <- lapply(unique(df$sel), \(x) {
    tmp <- df[df$sel == x, ]
    p <- ggplot(tmp, aes(cms_g, har_g, col = "#F4A582")) +
        geom_point(alpha = 0.2, size = 0.6) + 
        xlab(paste0("condition cms of ", x)) +
        ylab("condition cms of HVG + harmony") +
        facet_grid(t ~ s, labeller = \(.) label_both(.)) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        coord_fixed() + theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank()) + 
        ggtitle(x)
})

p2 <- lapply(unique(df$sel), \(x) {
    tmp <- df[df$sel == x, ]
    p <- ggplot(tmp, aes(cms_k, har_k)) +
        geom_point(alpha = 0.2, size = 0.6) + 
        xlab(paste0("cell type cms of ", x)) +
        ylab("cell type cms of HVG + harmony") +
        facet_grid(t ~ s, labeller = \(.) label_both(.)) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        coord_fixed() + theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank()) + 
        ggtitle(x)
})


tmp <- df[df$sel == df$sel[1], ]
p3 <- ggplot(tmp, aes(harmony_1, harmony_2, 
        col = cluster_re, shape = group_id)) +
        geom_point(alpha = 0.2, size = 0.6) + 
        facet_grid(t ~ s, labeller = \(.) label_both(.)) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        coord_fixed() + theme_bw(9) + theme(
            legend.position = "none",
            panel.grid = element_blank()) + 
        ggtitle("HVG + harmony")


plt <- append(p1, p2)

pdf(args[[2]], width = 10, height = 10, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
print(p3)
dev.off()