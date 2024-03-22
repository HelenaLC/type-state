#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-nCluster.pdf")

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
df <- lapply(res, select, 
  das, dat, sel, prop, nCells) |>
    do.call(what=rbind) |>
    distinct(das, dat, sel, nCells, prop) |>
    group_by(das, dat, sel) |>
    summarise_at("prop", median) |>
    mutate(sel=factor(sel, SEL))

aes <- list(
  scale_fill_gradientn(
    expression("median proportion"), 
    na.value="lightgrey", 
    colors=c("ivory", "pink", "red", "firebrick", "black")),
  labs(x="Feature selection", y="DS method"),
  coord_fixed(expand=FALSE),
  theme_minimal(6), 
  theme(
    plot.margin=margin(),
    panel.grid=element_blank(),
    legend.title=element_text(vjust=1),
    panel.border=element_rect(fill=NA),
    legend.key.width=unit(1, "lines"),
    legend.key.height=unit(0.5, "lines"),    
    plot.tag=element_text(size=9, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1)))


# plotting
gg <- ggplot(df) +
  geom_tile(
    aes(sel, das, fill=prop), 
    col="white", linewidth=0.1) +
  facet_grid(~dat) +
  aes

# saving
ggsave(args[[2]], gg, width=13, height=6, units="cm")