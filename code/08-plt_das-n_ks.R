#args <- list(list.files("outs/sim", "^das.*", full.names = TRUE), "plts/sim/das-n_ks.pdf")

# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
th <- 10
df <- lapply(res, select, t, s, 
  das, sel, i_nhood) |>
  do.call(what=rbind) |>
  distinct(t, s, das, sel, i_nhood) |>
  group_by(t, s, das, sel) |> 
  summarize(n=n(), .groups="drop") |>
  mutate(sel=factor(sel, c(DES, SEL))) |>
  mutate(n=case_when(n > th ~ th, TRUE ~ n)) 
j <- !(i <- df$sel %in% DES)
df_des <- df[i, ]
df_sel <- df[j, ]

# aesthetics
aes <- list(
  geom_tile(
    aes(t, s, fill=n), 
    col="white", linewidth=0.1),
  facet_grid(sel~das),
  scale_fill_gradientn("# clusters\n/nhoods.", 
    limits=c(0, th), breaks=c(0, th), 
    labels=\(.) ifelse(. == th, paste0(th, "+"), .),
    na.value="lightgrey", colors=c("ivory", "pink", "red", "firebrick", "black")),
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
