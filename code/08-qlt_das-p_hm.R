#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-p_hm.pdf")

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
df <- lapply(res, select, dat, 
  das, sel, i_nhood, p_nhood) |>
  do.call(what=rbind) |>
  distinct(dat, das, sel, 
    i_nhood, .keep_all=TRUE) |>
  group_by(dat, das, sel) |>
  summarise_at("p_nhood", mean) |>
  mutate(sel=factor(sel, SEL))

# aesthetics
ls <- range(df$p_nhood)
ls <- c(
  floor(ls[1]/(. <- 0.1))*., 
  ceiling(ls[2]/.)*.)

# plotting
gg <- ggplot(df) + facet_grid(~dat) +
  geom_tile(
    aes(sel, das, fill=p_nhood), 
    col="white", linewidth=0.1) +
  scale_fill_gradientn(
    "mean group prop.\nper cluster/nhood.",
    na.value="lightgrey", limits=ls, breaks=ls,
    colors=c("ivory", "pink", "red", "firebrick", "black")) +
  labs(x="Feature selection", y="DS method") +
  coord_fixed(expand=FALSE) +
  theme_minimal(6) + theme(
    plot.margin=margin(),
    panel.grid=element_blank(),
    legend.title=element_text(vjust=1),
    panel.border=element_rect(fill=NA),
    legend.key.size=unit(0.5, "lines"),    
    plot.tag=element_text(size=9, face="bold"),
    axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[2]], gg, width=7, height=3, units="cm")
