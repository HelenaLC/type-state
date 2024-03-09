# args <- list(list.files("outs/sim", "^das-", full.names=TRUE), "plts/sim/das-nCluster.pdf")

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
  t, s, b, sel, das, nSubpop) |>
  do.call(what=rbind) |>
  distinct(t, s, b, sel, das, nSubpop)

# separate ground truth selections
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

aes <- list(geom_tile(col="white"),
    facet_grid(t~s, scales="free"),
    labs(x="feature selection", y="DS methods"),
    geom_text(aes(label = nSubpop), size = 2.5),
    theme(
      plot.margin=margin(),
      panel.grid=element_blank(),
      legend.title=element_text(vjust=1),
      axis.text.x=element_text(angle=45, hjust=1, vjust=1)),
      scale_fill_distiller(NULL, palette="YlOrRd", direction = 1))

p1 <- ggplot(df_des, aes(sel, das, fill=as.numeric(nSubpop))) + aes 
p2 <- ggplot(df_sel, aes(sel, das, fill=as.numeric(nSubpop))) + aes 
  
pdf(args[[2]], onefile=TRUE, width=15, height=16)
for (p in list(p1,p2)) print(p); dev.off()