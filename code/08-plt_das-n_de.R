#args <- list(list.files("outs/sim", "^das-", full.names=TRUE), "plts/sim/das-n_de.pdf")

# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggrastr)
  library(ggplot2)
  library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
th <- 400
df <- lapply(res, select, t, s, 
    sel, das, gene_id, p_adj) |>
    do.call(what=rbind) |>
    group_by(t, s, sel, das, gene_id) |>
    summarise_at("p_adj", min, na.rm=TRUE) |>
    group_by(t, s, sel, das) |>
    summarize(n_de=sum(p_adj < 0.05), .groups="drop") |>
    mutate(n_de=case_when(n_de > th ~ th, TRUE ~ n_de)) |>
    complete(t, s, sel, das, fill=list(n_de=NA)) |>
    mutate(das=factor(gsub("^DS_", "", das), DAS)) |>
    mutate(sel=factor(sel, c(DES, SEL)))

# split by ground truth-/
# method-based selection
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# plotting
aes <- list(
    geom_tile(
        aes(t, s, fill=n_de), 
        col="white", linewidth=0.1),
    facet_grid(sel~das),
    scale_fill_gradientn(
        labels=\(.) ifelse(. == th, paste0(th, "+"), .),
        "# DEGs (any adj.\np-value < 0.05)", trans="sqrt",
        limits=c(0, th), breaks=c(0, 50, 100, seq(200, th, 200)), 
        na.value="lightgrey", colors=c("ivory", "gold", "red", "navy")),
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

gg <- ggplot(df_des) + ggplot(df_sel) + 
  plot_layout(nrow=1, guides="collect") & 
  plot_annotation(tag_levels="a") &
  aes & theme(
    plot.margin=margin(),
    legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=12.5, units="cm")
