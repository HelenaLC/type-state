# args <- list(
#     list.files("outs", "^das.*", full.names = TRUE),
#     "plts/das-p_val_density.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(stringr)
    library(dplyr)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(.) .[c("sel","t", "s", "das", "p_val", "p_adj")])

typ <- vapply(res, \(.) .$das[1], character(1))
ds <- !(da <- grepl("^DA", typ))
da <- do.call(rbind, res[da])
ds <- do.call(rbind, res[ds])
all <- do.call(rbind, res)
all <- all %>% 
    mutate(
        DAS = ifelse(str_detect(das, "DS"), "DS", "DA"),
        clustering = ifelse(str_detect(das, "limma|edgeR"), "Yes", "No"))

p1 <- lapply(unique(all$sel), \(s){
    df <- all[all$sel == s,]
    gg <- ggplot(df, aes(p_val, after_stat(ndensity), col = das, linetype = DAS)) +
        facet_grid(s ~ t, labeller = \(.) label_both(.)) +
        geom_density(key_glyph = "smooth") +
        #guides(col = guide_legend(override.aes = list(size = 3)),
        #    linetype =  guide_legend(override.aes = list(linetype = 1))) +
        scale_x_continuous("p-value", breaks = seq(0, 1, 0.5), expand = c(0, 0.05)) +
        scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
        ggtitle(s) + 
        theme_bw(9) + theme(
            panel.grid = element_blank(),
            axis.line = element_line(size = 0.25),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.key.size = unit(0.5, "line"),
            legend.text = element_text(size = 7), 
            axis.text.x = element_text(angle = 45, hjust = 1)) 
})



# out <- lapply(args[[1]], \(x) { if (str_detect(x,"dat")) readRDS(x) })
# out  <- out[!vapply(out, is.null, logical(1))]
# out <- lapply(out, \(.) .[c("sel", "dat", "das", "p_val", "p_adj")])
# out <- do.call(rbind, out)
# out <- out %>% 
#     mutate(
#         DAS = ifelse(str_detect(das, "DS"), "DS", "DA"),
#         clustering = ifelse(str_detect(das, "limma|edgeR"), "Yes", "No"))
# 
# p2 <- lapply(unique(out$dat), \(s){
#     fd <- out[out$dat == s,]
#     p <- ggplot(fd, aes(p_val, after_stat(ndensity), col = das, linetype = DAS)) +
#         geom_density(key_glyph = "smooth") +
#         facet_wrap(~ sel, ncol = 3, scales = "free") + 
#         scale_x_continuous("p-value", breaks = seq(0, 1, 0.5), expand = c(0, 0.05)) +
#         scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
#         ggtitle(s) + 
#         theme_bw(9) + theme(
#             panel.grid = element_blank(),
#             axis.line = element_line(size = 0.25),
#             plot.title = element_text(hjust = 0.5, size = 15),
#             legend.key.size = unit(0.5, "line"),
#             legend.text = element_text(size = 7), 
#             axis.text.x = element_text(angle = 45, hjust = 1)) 
# })

#plt <- append(p1,p2)
plt <- p1


pdf(args[[2]], width = 12, height = 8, onefile = TRUE)
for (i in seq_along(plt)) {
    print(plt[[i]])
}
dev.off()



