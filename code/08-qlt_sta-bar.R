suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)
# wrangling
.f <- \(df) df %>%
  group_by(dat, sel, sta) %>%
  summarise_at("sta_val", mean) 
df <- .f(df)

gg <- ggplot(df, aes(x = sel, y = sta_val, color = sel)) +
  geom_bar(stat="identity", fill="white", position=position_dodge())+
  xlab("Evaluation metrics") +
  ylab("Metric values") +
  labs(color="Feature selection") + 
  facet_grid2(dat ~ sta, scales = "free", independent = "y") +
  theme_bw() +
  scale_color_brewer(palette = "Paired") & theme(
    legend.position="bottom",
    legend.justification=c(0.5, 1),
    legend.box.spacing=unit(0, "pt"),
    panel.grid.minor=element_blank(), 
    panel.border=element_rect(fill=NA),
    plot.tag=element_text(size=9, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1))


ggsave(args[[2]], gg, width=30, height=15, units="cm")

