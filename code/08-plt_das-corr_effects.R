# args <- list(
#     list.files("outs/sim", "^das-", full.names=TRUE),
#     "plts/sim/das-corr_effects.pdf")

# dependencies
suppressPackageStartupMessages({
    library(poolr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- lapply(res, select, t, s,
    sel, das, p_adj, gene_id) |>
    do.call(what=rbind) |>
    filter(!is.na(p_adj)) |>
    group_by(t, s, sel, das, gene_id) |>
    summarize_at("p_adj", ~fisher(.x)$p)

# simulations
rd <- list.files("data/sim/01-pro", "-rd\\.rds", full.names=TRUE)
rd <- do.call(rbind, lapply(rd, readRDS))

de <- grep("^GroupDE", names(rd))
ds <- grep("^ConditionDE", names(rd))

rd$de <- sapply(de, \(i) {
    j <- setdiff(de, i)
    lfc <- log2(rd[[i]]/rd[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

rd$ds <- sapply(ds, \(i) {
    j <- setdiff(ds, i)
    lfc <- log2(rd[[i]]/rd[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

# merging
rd <- rd[c("t", "s", "gene_id", "de", "ds")]
df <- merge(x=df, y=rd, by=c("t", "s", "gene_id"))

# correlation
dfs <- split(df, df[. <- c("t", "s", "sel", "das")])
fd <- lapply(dfs, \(df) {
    t <- cor(df$p_adj, df$de, method="spearman")
    s <- cor(df$p_adj, df$ds, method="spearman")
    data.frame(
        row.names=NULL, df[1, .], 
        cor=c("cor_t", "cor_s"), 
        cor_val=c(t, s)) 
}) |> 
    do.call(what=rbind) |> na.omit() |> 
    mutate(das=gsub("^DS_", "", das))

# split selections
fd <- mutate(fd, sel=factor(sel, c(DES, SEL)))

. <- filter(fd, cor == "cor_t")
j <- !(i <- .$sel %in% DES)
df_t <- .[i, ]; fd_t <- .[j, ]
rng_t <- range(.$cor_val, na.rm=TRUE)

. <- filter(fd, cor == "cor_s")
j <- !(i <- .$sel %in% DES)
df_s <- .[i, ]; fd_s <- .[j, ]
rng_s <- range(.$cor_val, na.rm=TRUE)

# aesthetics
rng_t <- c(floor(rng_t[1]*10)/10, ceiling(rng_t[2]*10)/10)
rng_s <- c(floor(rng_s[1]*10)/10, ceiling(rng_s[2]*10)/10)

pal_t <- scale_fill_gradient2(limits=rng_t, #breaks=c(rng_t, 0),
    na.value="lightgrey", low="blue", mid="ivory", high="red")
pal_s <- scale_fill_gradient2(limits=rng_s, #breaks=c(rng_s, 0),
    na.value="lightgrey", low="blue", mid="ivory", high="red")

aes <- list(
    facet_grid(das ~ sel),
    coord_fixed(expand=FALSE),
    scale_x_continuous("type effect", n.breaks=2),
    scale_y_continuous("state effect", n.breaks=2),
    geom_tile(col="white", linewidth=0.1, aes(t, s, fill=cor_val)),
    theme_minimal(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(vjust=1),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines")))

# plotting
lab_t <- expression("Cor("*X^2*", logFC"[type]*")")
lab_s <- expression("Cor("*X^2*", logFC"[state]*")")
p1 <- 
    ggplot(df_t) + ggplot(fd_t) + 
    plot_layout(ncol=1, guides="collect") &
    plot_annotation(tag_levels="a") &
    pal_t & labs(fill=lab_t) & 
    aes & theme(
        plot.margin=margin(), 
        legend.position="bottom",
        plot.tag=element_text(size=9, face="bold"))
p2 <- 
    ggplot(df_s) + ggplot(fd_s) +
    plot_layout(ncol=1, guides="collect") &
    plot_annotation(tag_levels="a") &
    pal_s & labs(fill=lab_s) & 
    aes & theme(
        plot.margin=margin(), 
        legend.position="bottom",
        plot.tag=element_text(size=9, face="bold"))

# saving
pdf(args[[2]], onefile=TRUE, width=12/2.54, height=14/2.54)
for (p in list(p1, p2)) print(p); dev.off()



