#args <- list(list.files("outs/sim", "^das-", full.names=TRUE), "plts/sim/das-p_adj.pdf")

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
df <- lapply(res, select, t, s, 
  das, sel, gene_id, p_adj) |>
  do.call(what=rbind) |>
  group_by(t, s, das, sel, gene_id) |>
  summarise_at("p_adj", min, na.rm=TRUE) |>
  mutate(
    t=factor(t, seq(0, 1, 0.2)),
    s=factor(s, seq(0, 1, 0.2)),
    sel=factor(sel, c(DES, SEL)), 
    das=factor(gsub("^DS_", "", das), DAS))
foo <- c("random", "DEnotDS", "DS")
fd <- df |>
  filter(sel %in% foo) |>
  mutate(sel=factor(sel, foo))

# plotting
gg <- ggplot(fd, 
  aes(das, p_adj, col=t)) +
  facet_grid(s~sel, labeller=labeller(
    sel=label_value, s=label_both)) +
  geom_boxplot(
    key_glyph="point", width=0.75,
    outlier.size=0.5, linewidth=0.25,
    outlier.shape=16, outlier.stroke=0) +
  geom_hline(yintercept=0.05, lty=2, linewidth=0.1) +
  scale_color_brewer("type\neffect", palette="Blues") +
  guides(color=guide_legend(override.aes=list(size=2))) +
  labs(x="DSA method", y="min. adj. p-value (across clusters/nhoods.)") +
  scale_y_continuous(expand=c(0, 0.1), limits=c(0, 1), n.breaks=2) +
  scale_x_discrete(expand=c(0, 0.5)) +
  theme_minimal(6) + theme(
    plot.margin=margin(),
    panel.grid=element_blank(),
    panel.border=element_rect(fill=NA),
    legend.key.size=unit(0.5, "lines"),
    legend.title=element_text(vjust=1.5),
    plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=10, height=10, units="cm")
