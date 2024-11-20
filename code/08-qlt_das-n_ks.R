#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-n_ks.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
th <- 100
fd <- lapply(res, select,
  dat, sel, das, i_nhood) |>
  do.call(what=rbind) |>
  distinct(dat, sel, das, i_nhood) |>
  group_by(dat, sel, das) |> 
  summarize(n=n(), .groups="drop") |>
  mutate(sel=factor(sel, SEL))
df <- mutate(fd, n=case_when(n > th ~ th, TRUE ~ n))

# plotting
gg <- ggplot(df, aes(das, sel, fill=n)) + 
    geom_tile(col="white") + facet_grid(~dat) +
    geom_text(data=fd, aes(label=n), size=1.5) +
    scale_fill_gradientn("# clusters\n/nhoods.", 
        limits=c(0, th), labels=\(.) ifelse(. == th, paste0(th, "+"), .),
        colors=c("ivory", "pink", "red", "firebrick"), na.value="lightgrey") +
    labs(x="feature selection", y="DSA method") +
    coord_equal(expand=FALSE) +
    theme_minimal(6) + theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        legend.title=element_text(vjust=1),
        legend.key.size=unit(0.5, "lines"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[2]], gg, width=4, height=4, units="cm")
