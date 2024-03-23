# args <- list(
#     list.files("outs/sim", "^das-", full.names=TRUE),
#     "plts/sim/das-F1.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(caret)
    library(patchwork)
    library(matrixStats)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- lapply(res, select, das, sel, p_adj, gene_id, t, s,
    starts_with("GroupDE.Group"), 
    starts_with("ConditionDE.Condition")) |>
    do.call(what=rbind) |>
    group_by(das, sel, gene_id, t, s,
        across(starts_with("GroupDE.Group")),
        across(starts_with("ConditionDE.Condition")))  |>
    summarize_at("p_adj", min) |>
    mutate(test_ds = p_adj < 0.05) |>
    ungroup()

de <- grep("^GroupDE", names(df))
ds <- grep("^ConditionDE", names(df))
DE <- rowAnys(df[de] != 1)
DS <- rowAnys(df[ds] != 1)
df$DS <- DS
df$DSnotDE <- DS & !DE


DE <- sapply(de, \(i) {
    j <- setdiff(de, i)
    lfc <- log2(df[j]/df[[i]])
    rowMeans(lfc) }) |> rowMeans()

DS <- sapply(ds, \(i) {
    j <- setdiff(ds, i)
    lfc <- log2(df[j]/df[[i]])
    rowMeans(lfc) }) |> rowMeans()

df$DSgtDE <- DS > DE
ls <- c("TRUE", "FALSE")
df <- df %>%
    group_by(das, sel, t, s) %>%
    mutate(F1_DS = confusionMatrix(factor(test_ds, levels=ls), 
        factor(DS,levels=ls), positive = "TRUE")$byClass["F1"]) %>%
    mutate(F1_DSnotDE = confusionMatrix(factor(test_ds, levels=ls), 
        factor(DSnotDE,levels=ls), positive = "TRUE")$byClass["F1"]) %>%
    mutate(F1_DSgtDE = confusionMatrix(factor(test_ds, levels=ls), 
        factor(DSgtDE,levels=ls), positive = "TRUE")$byClass["F1"])

aes <- list(
    geom_tile_rast( 
        col="white", linewidth=0.1),
    facet_grid(sel~das),
    scale_fill_gradientn(
        "F1",
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

df <- mutate(df, sel=factor(sel, c(DES, SEL)))
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# plotting
p1 <- ggplot(df_des, aes(t, s, fill=F1_DS)) + 
    ggplot(df_sel, aes(t, s, fill=F1_DS)) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(title="F1_DS", tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

p2 <- ggplot(df_des, aes(t, s, fill=F1_DSnotDE)) + 
    ggplot(df_sel, aes(t, s, fill=F1_DSnotDE)) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(title="F1_DSnotDE", tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

p3 <- ggplot(df_des, aes(t, s, fill=F1_DSgtDE)) + 
    ggplot(df_sel, aes(t, s, fill=F1_DSgtDE)) + 
    plot_layout(nrow=1, guides="collect") & 
    plot_annotation(title="F1_DSgtDE", tag_levels="a") &
    aes & theme(
        plot.margin=margin(),
        legend.position="bottom")

# saving
pdf(args[[2]], onefile=TRUE, width=6, height=6.5)
for (p in list(p1,p2,p3)) print(p); dev.off()


