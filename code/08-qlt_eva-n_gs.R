#args <- list(list.files("outs/dat", "^eva-", full.names=TRUE), "plts/dat/eva-line.pdf")

# dependencies
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggh4x)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- do.call(rbind, res) |>
  distinct(dat, num, sta, sel, sta_val, .keep_all=TRUE) %>%
  group_by(dat, num, sta, sel, nGenes) %>% 
  summarize(sta_val=mean(sta_val)) %>% 
  mutate(sel=factor(sel, SEL)) 

# plotting
pal <- hcl.colors(nlevels(df$sel), "Spectral")
gg <- ggplot(df, aes(nGenes, sta_val, col=sel, group=sel)) + 
  facet_wrap(~sta, ncol=4, scales="free") +
  geom_vline(xintercept=2e3, lty=2, linewidth=0.2) +
  geom_point(alpha=0.8, size=0.5) +
  geom_line(alpha=0.8, linewidth=0.5, show.legend=FALSE) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=2))) +
  labs(x="# selected features", y="statistic value") + 
  scale_color_manual("selection", values=pal) +
  theme_bw(6) + theme(
    legend.position="bottom",
    panel.grid=element_blank(), 
    legend.justification=c(0.5, 1),
    legend.box.spacing=unit(0, "pt"),
    panel.border=element_rect(fill=NA),
    legend.key.size=unit(0.5, "lines"),
    plot.tag=element_text(size=9, face="bold"))

# saving
ggsave(args[[2]], gg, width=15, height=6, units="cm")
