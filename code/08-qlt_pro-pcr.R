args <- list(list.files("data/dat/02-rep", "Kang18.*-cd", full.names=TRUE), "plts/dat/pro-pcr.pdf")

# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(SingleCellExperiment)
})

# loading
dfs <- lapply(args[[1]], readRDS)

# analysis
df <- lapply(dfs, \(df) {
    y <- df[, grep("^PC[0-9]+$", names(df))]
    z <- lapply(c("group_id", "cluster_id", "cluster_re"), \(x) {
        z <- summary(lm(as.matrix(y) ~ df[[x]]))
        r2 <- sapply(z, \(.) .$adj.r.squared)
        data.frame(x, pc=seq_along(r2), r2, sel=df$sel[1])
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# wrangling
df <- df |>
    mutate(sel=factor(sel, SEL))

# plotting
gg <- ggplot(df, aes(pc, r2, col=x)) +
    coord_cartesian(xlim=c(1, 16)) +
    facet_grid(~sel) + geom_point(size=0.6) +
    geom_line(linewidth=0.3, show.legend=FALSE) +
    scale_x_continuous(breaks=c(1, seq(5, 100, 5))) +
    scale_y_continuous(limits=c(0, 1), n.breaks=3) +
    guides(col=guide_legend(override.aes=list(size=2))) +
    scale_color_manual("predictor", values=c("cyan", "gold", "magenta")) +
    labs(x="principal component", y="coeff. of determination\nfrom linear regression") + 
    theme_minimal(6) + theme(
        aspect.ratio=1,
        panel.grid=element_blank(), 
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.5, "lines"))

# saving
ggsave(args[[2]], gg, width=16, height=3.5, units="cm")
