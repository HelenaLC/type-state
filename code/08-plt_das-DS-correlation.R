# args <- list(
#     list.files("outs", "^das.*", full.names = TRUE),
#     "plts/das-p_val_density.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(stringr)
    library(dplyr)
    library(poolr)
    library(reshape2)
    library(GGally)
    library(ggrastr)
    library(tidyr)
    library(ggpubr)
    library("viridis")
})

# read DS results
res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim") & str_detect(x,"DS")) 
    readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]
res <- lapply(res, \(.) .[c("sel","t", "s", "das", "p_val", 
    "p_adj", "gene_id")])
tbl <- lapply(res, \(x) {
    out <- vapply(unique(x$gene_id), \(g) {
        p <- x[x$gene_id == g, "p_adj"]
        #fisher(p)$statistic/attr(fisher(p)$statistic,"df")
        -log(fisher(p)$p)
    }, numeric(1))
    data.frame(row.names = NULL, gene_id = names(out), fisher = out,
        t = x$t[seq_len(length(out))], s = x$s[seq_len(length(out))], 
        sel = x$sel[seq_len(length(out))], das = x$das[seq_len(length(out))])
})
ds <- do.call(rbind, tbl)

# read ground truth
files <- list.files("data/01-fil", "*rd.rds", full.names = TRUE)
rd <- lapply(files, \(x) readRDS(x))
df <- do.call(rbind, rd)

## define type effects
gde <- df[grep("GroupDE", names(df))]
mg <- sapply(seq_len(ncol(gde)), \(i){
    not_i <- setdiff(seq_len(ncol(gde)), i)
    mgk <- sapply(not_i, \(j) {
        log(gde[,i]/gde[,j], base = 2)
    })
    rowMeans(mgk)
})

df$mg <- apply(abs(mg), 1, mean)

## define state effects
cde <- df[grep("ConditionDE", names(df))]
cg <- sapply(seq_len(ncol(cde)), \(i){
    not_i <- setdiff(seq_len(ncol(cde)), i)
    cgk <- sapply(not_i, \(j) {
        log(cde[,i]/cde[,j], base = 2)
    })
    rowMeans(cgk)
})
df$cg <- apply(abs(cg), 1, mean)
df <- df[,c("t", "s", "gene_id", "cg", "mg")]

all <- merge(x = df, y = ds, by = c("t", "s", "gene_id"))

# correlation
lst <- split(all, list(all$das, all$t, all$s, all$sel))
scr <- lapply(lst, \(x) {
    cor_s <- cor(x$fisher, x$cg, method = "spearman")
    cor_t <- cor(x$fisher, x$mg, method = "spearman")
    rbind(data.frame(das = x$das[1], t = x$t[1], 
        s = x$s[1], sel = x$sel[1], cor_val = cor_s, cor = "cor_state"),
        data.frame(das = x$das[1], t = x$t[1], 
            s = x$s[1], sel = x$sel[1], cor_val = cor_t, cor = "cor_type"))
}) 
#scr <- scr[!is.na(scr)]
res <- do.call(rbind, scr)
res <- na.omit(res)
res_s <- res[res$cor == "cor_state", ]
res_t <- res[res$cor == "cor_type", ]

p1 <- ggplot(res_s) +
    facet_grid(das ~ sel) +
    geom_tile(col = "white", aes(t, s, fill = cor_val)) +
    scale_x_continuous() +
    scale_y_continuous() +
    labs(x = "type effect", y = "state effect") +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        limits=c(-1, 1), n.breaks=3, direction=-1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA)) +
    ggtitle("Spearman's correlation with state effect")

p2 <- ggplot(res_t) +
    facet_grid(das ~ sel) +
    geom_tile(col = "white", aes(t, s, fill = cor_val)) +
    scale_x_continuous() +
    scale_y_continuous() +
    labs(x = "type effect", y = "state effect") +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        limits=c(-1, 1), n.breaks=3, direction=-1) +
    coord_fixed(expand = FALSE) + 
    theme_minimal(9) + theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA)) +
    ggtitle("Spearman's correlation with type effect")

# correlation between methods
ndf <- dcast(data = all,formula = ...~das,
    fun.aggregate = sum, 
    value.var = "fisher")

# plt <- lapply(unique(ndf$sel), \(x) {
#     df <- ndf[ndf$sel == x,]
#     p <- ggpairs(df, columns = c("DS_lemur", "DS_edgeR", "DS_miloDE")) +
#         ggtitle(x)
# }) 

# scatter plot between methods
fd <- all |>
    group_by(t, s, sel, das) |>
    # mutate(
    #     fisher=fisher-min(fisher),
    #     fisher=fisher/max(fisher)) |>
    ungroup() |>
    mutate(
        across(all_of(c("t", "s")), 
            ~factor(., sort(unique(.))))) |>
    pivot_wider(
        names_from="das",
        values_from="fisher",
        id_cols=setdiff(names(all), c("das", "fisher"))) |>
    pivot_longer(
        names_to="effect",
        cols=c("cg", "mg")) |>
    mutate(effect=factor(effect, 
        levels=c("cg", "mg"), 
        labels=c("state", "type"))) |>
    ungroup() |>
    arrange(value)

pa <- ggplot(fd, aes(DS_edgeR, DS_miloDE, col=value)) 
pb <- ggplot(fd, aes(DS_edgeR, DS_lemur, col=value))
#pc <- ggplot(fd, aes(DS_miloDE, DS_lemur, col=value))

p3 <- wrap_plots(pa, pb, ncol=1) +
    plot_layout(guides="collect") &
    facet_grid(effect~sel) &
    geom_point_rast(shape=16, alpha=0.1, size=0.01) &
    scale_color_viridis(option = "C") &
    scale_x_continuous(n.breaks=3) &
    scale_y_continuous(n.breaks=3) &
    stat_cor(method="pearson", aes(label = ..r.label..)) &
    theme_linedraw(6) & theme(
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        legend.key.size=unit(0.5, "lines"))


pdf(args[[2]], width = 15, height = 7, onefile = TRUE)
# for (i in seq_along(plt)) {
#     print(plt[[i]])
# }
print(p1)
print(p2)
print(p3)
dev.off()
