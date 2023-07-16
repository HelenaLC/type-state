suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
    library(stringr)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    mutate(sco_val = case_when(
        !(sco %in% c("type_entropy","random")) ~ log10(sco_val),
        .default = sco_val)) 

# DE genes only (according to 'Splatter')
de <- df[grep("GroupDE", names(df))]
fd <- df[!rowAlls(as.matrix(de) == 1), ]

gg <- list(
    geom_violin(trim = FALSE, fill = "gray"),
    geom_boxplot(width = 0.1),
    guides(color = guide_legend(override.aes = list(size = 2))))

labs <- \(.) label_value(., multi_line = FALSE)
#facet_grid( ~ s, labeller = \(.) label_both(.), scales = "free")
# facet_wrap(t ~ sco, labeller = labs, scales = "free") 
fac_s <- facet_wrap(s ~ sco, labeller = labs, scales = "free") 
fac_t <- facet_wrap(t ~ sco, labeller = labs, scales = "free") 

pal_s <- scale_color_brewer(palette = "Reds", "state\neffect", limits = seq(0, 1, 0.2))
pal_t <- scale_color_brewer(palette = "Blues", "type\neffect", limits = seq(0, 1, 0.2))

p1 <- ggplot(df, aes(x = factor(t), y = sco_val, col = factor(t))) + gg + fac_s + pal_t
p2 <- ggplot(fd, aes(x = factor(t), y = sco_val, col = factor(t))) + gg + fac_s + pal_t

p3 <- ggplot(df, aes(x = factor(s), y = sco_val, col = factor(s))) + gg + fac_t + pal_s
p4 <- ggplot(fd, aes(x = factor(s), y = sco_val, col = factor(s))) + gg + fac_t + pal_s

thm <- theme_bw(6) + theme(
    panel.grid = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing = unit(2, unit = "mm"),
    plot.margin = margin(0, unit = "mm"),
    legend.key.size = unit(0.5, "lines"),
    strip.text = element_text(color = "black"),
    plot.tag = element_text(size = 9, face = "bold"),
    strip.background = element_rect(color = NA, fill = "white"))

f1 <- (p1 / p3) & plot_annotation(tag_levels = "a") & thm
f2 <- (p2 / p4) & plot_annotation(tag_levels = "a") & thm

dat <- lapply(args[[1]], \(x) { if (str_detect(x,"dat")) readRDS(x) })
dat <- dat[!vapply(dat, is.null, logical(1))]

rdf <- do.call(rbind, dat) %>%
    mutate(sco_val = case_when(
        !(sco %in% c("type_entropy","random")) ~ log10(sco_val),
        .default = sco_val)) 

plt <- lapply(unique(rdf$dat), \(x) {
    tmp <- rdf[rdf$dat == x,]
    p <- ggplot(tmp, aes(x = sco, y = sco_val)) + 
        gg + ggtitle(x) 
})
names(plt) <- unique(rdf$dat)



pdf(args[[2]], width = 15/1.5, height = 18/1.5, onefile = TRUE)
print(f1) 
print(f2)
for (x in names(plt)) print(plt[[x]]) 
dev.off()

