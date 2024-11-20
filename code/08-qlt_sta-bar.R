args <- list(list.files("outs/dat", "^sta-", full.names=TRUE), "plts/dat/sta-bar.pdf")

# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidytext)
  library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
df <- do.call(rbind, res) 
.f <- \(df) df %>%
    group_by(dat, sel, sta) %>%
    summarise_at("sta_val", mean) %>% 
    mutate(sel=factor(sel, SEL))
df <- .f(df)

# plotting
gg <- ggplot(df, aes(y=sta_val, fill=sel,
    x=reorder_within(sel, sta_val, sta))) +
    facet_wrap(~sta, scales="free", nrow=1) +
    geom_bar(width=0.8, stat="identity", key_glyph="point") +
    scale_fill_brewer("feature\nselection\nmethod", palette="Set2") + 
    guides(fill=guide_legend(override.aes=list(shape=21, size=2, stroke=0))) + 
    scale_y_continuous(NULL, expand=expansion(0.02, 0)) +
    scale_x_reordered(NULL, expand=expansion(0, 0.6)) +
    theme_minimal(6) + theme(
        aspect.ratio=1.5,
        panel.grid=element_blank(), 
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.5, "lines"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[2]], gg, width=15, height=3.5, units="cm")
