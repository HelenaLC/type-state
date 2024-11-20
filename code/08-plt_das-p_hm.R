#args <- list(list.files("outs/sim", "^das.*", full.names = TRUE), "plts/sim/das-p_hm.pdf")

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
  das, sel, i_nhood, p_nhood) |>
  do.call(what=rbind) |>
  distinct(t, s, das, sel, 
    i_nhood, .keep_all=TRUE) |>
  group_by(t, s, das, sel) |>
  summarize_at("p_nhood", mean) |>
  mutate(sel=factor(sel, c(DES, SEL)))
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# aesthetics
ls <- range(df$p_nhood, na.rm=TRUE)
ls <- c(floor(ls[1]/(. <- 0.1))*., ceiling(ls[2]/.)*.)
aes <- list(
  geom_tile(
    aes(t, s, fill=p_nhood), 
    col="white", linewidth=0.1),
  facet_grid(sel~das),
  scale_fill_gradientn(
    "mean group prop.\nper cluster/nhood.",
    limits=ls, breaks=ls, na.value="lightgrey", 
    colors=c("ivory", "pink", "red", "firebrick", "black")),
  labs(x="type effect", y="state effect"),
  scale_x_continuous(breaks=c(0, 1)),
  scale_y_continuous(breaks=c(0, 1)),
  coord_fixed(expand=FALSE),
  theme_minimal(6), theme(
    plot.margin=margin(),
    panel.grid=element_blank(),
    axis.title=element_text(hjust=0),
    panel.border=element_rect(fill=NA),
    legend.key.size=unit(0.5, "lines"),
    legend.title=element_text(vjust=1.5),
    plot.tag=element_text(size=9, face="bold")))

# plotting
gg <- ggplot(df_des) + ggplot(df_sel) + 
  plot_layout(nrow=1, guides="collect") & 
  plot_annotation(tag_levels="a") &
  aes & theme(
    plot.margin=margin(),
    legend.position="bottom")

# saving
ggsave(args[[2]], gg, width=12, height=12.5, units="cm")
