# args <- list(
#     list.files("outs/sim", "^das-", full.names=TRUE),
#     "plts/sim/das-F1.pdf")

# dependencies
suppressPackageStartupMessages({
    library(caret)
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
th <- 0.05
de <- "GroupDE"
ds <- "ConditionDE"

df <- lapply(res, select, t, s,
  das, sel, gene_id, p_adj, 
  starts_with(c(de, ds))) |>
  do.call(what=rbind) |>
  group_by(t, s, das, sel, gene_id) |>
  summarise(
    res=any(p_adj < th), .groups="drop",
    across(starts_with(c(de, ds)), unique))

de_idx <- grep("^GroupDE", names(df))
ds_idx <- grep("^ConditionDE", names(df))
is_de <- rowAnys(df[de_idx] != 1)
is_ds <- rowAnys(df[ds_idx] != 1)
df$DSnotDE <- (df$DS <- is_ds) & !is_de

de_lfc <- sapply(de_idx, \(i) {
    j <- setdiff(de_idx, i)
    lfc <- log2(df[j]/df[[i]])
    rowMeans(lfc) }) |> rowMeans()

ds_lfc <- sapply(ds_idx, \(i) {
    j <- setdiff(ds_idx, i)
    lfc <- log2(df[j]/df[[i]])
    rowMeans(lfc) }) |> rowMeans()

df$DSgtDE <- ds_lfc > de_lfc

ls <- c("TRUE", "FALSE")
df <- df |>
    mutate(das=gsub("^DS_", "", das)) |>
    mutate(das=factor(das, DAS)) |>
    group_by(das, sel, t, s) |>
    mutate(
        F1_DS=confusionMatrix(
            factor(res, levels=ls), 
            factor(DS, levels=ls), 
            positive="TRUE")$byClass["F1"],
        F1_DSgtDE=confusionMatrix(
            factor(res, levels=ls), 
            factor(DSgtDE,levels=ls), 
            positive="TRUE")$byClass["F1"],
        F1_DSnotDE=confusionMatrix(
            factor(res, levels=ls), 
            factor(DSnotDE, levels=ls),
            positive="TRUE")$byClass["F1"])

# split by ground truth- & method-based selection
df <- mutate(df, sel=factor(sel, c(DES, SEL)))
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# plotting
aes <- list(
    geom_tile_rast( 
        col="white", linewidth=0.1),
    facet_grid(sel~das),
    scale_fill_gradientn(
        colors=c("ivory", "gold", "red", "navy"),
        na.value="lightgrey", limits=c(0, 1), n.breaks=3),
    labs(x="type effect", y="state effect"),
    scale_x_continuous(breaks=c(0, 1)),
    scale_y_continuous(breaks=c(0, 1)),
    coord_fixed(expand=FALSE),
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.title=element_text(vjust=1.5),
        plot.tag=element_text(size=9, face="bold")))

fs <- grep("F1", names(df), value=TRUE)
lb <- sprintf("(any adj.\np-value < %s)", th)
ps <- lapply(fs, \(.) {
    ggplot(df_des, aes(t, s, fill=.data[[.]])) + 
    ggplot(df_sel, aes(t, s, fill=.data[[.]])) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(tag_levels="a") & 
    labs(fill=paste(., lb)) & 
        aes & theme(
        plot.margin=margin(), 
        legend.position="bottom")
})

# saving
pdf(args[[2]], onefile=TRUE, width=12/2.54, height=12.5/2.54)
for (p in ps) print(p); dev.off()
