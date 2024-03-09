#args <- list(list.files("outs/sim", "^das.*", full.names = TRUE), "plts/das-p_val_density.pdf")

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
  das, t, s, sel, prop) |>
  do.call(what=rbind) |>
  distinct(das, t, s, sel, prop)

# plotting
ps <- lapply(split(df, df$sel), \(fd)
  ggplot(fd, aes(das, prop, col=das)) + 
  geom_violin_rast(position = position_dodge(width = 0.8), trim = FALSE) + 
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position = "none") +
  ggh4x::facet_grid2(t ~ s, scales = "free", independent = "y",
    labeller=\(.) label_both(.)) +
  scale_color_brewer(palette = "Paired") +
  ggtitle(fd$sel[1]))

# saving
pdf(args[[2]], onefile=TRUE, width=12, height=10)
for (p in ps) print(p); dev.off()