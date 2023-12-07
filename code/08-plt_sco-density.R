#args <- list(list.files("outs/sim", "^sco-", full.names=TRUE), "plts/sim/sco-density.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

# loading
res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

# wrangling
ex <- c("random", "HVG")
sco_ord <- c(ex, "type_Fstat", "type_PVE", "state_PVE", "state_edgeR")
df <- bind_rows(res) |>
    #mutate(sco_val=case_when(!sco %in% ex ~ log10(sco_val), TRUE ~ sco_val)) |>
    mutate(sco=factor(sco, sco_ord))

# subsetting
de <- grep("GroupDE", names(df))
ds <- grep("ConditionDE", names(df))
fd <- df[rowAnys(df[de] != 1) & !rowAnys(df[ds] != 1), ]

# aesthetics
aes <- list(
    scale_x_continuous(NULL, n.breaks=3),
    scale_y_continuous("scaled density", n.breaks=2),
    guides(color=guide_legend(override.aes=list(size=1))),
    geom_density(linewidth=0.4, key_glyph="point"),
    theme_bw(6), theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        panel.spacing=unit(2, unit="mm"),
        legend.key.size=unit(0.25, "lines"),
        strip.text=element_text(color="black"),
        plot.tag=element_text(size=9, face="bold"),
        strip.background=element_rect(color=NA, fill="white")))

.labs <- \(.) label_both(., multi_line=FALSE)
fac_s <- facet_grid(s ~ sco, labeller=.labs, scales="free") 
fac_t <- facet_grid(t ~ sco, labeller=.labs, scales="free") 
pal_s <- scale_color_brewer(palette="Reds", "state\neffect", limits=seq(0, 1, 0.2))
pal_t <- scale_color_brewer(palette="Blues", "type\neffect", limits=seq(0, 1, 0.2))

# plotting
p1 <- ggplot(df, aes(x=sco_val, y=..scaled.., col=factor(t))) + fac_s + pal_t
p2 <- ggplot(df, aes(x=sco_val, y=..scaled.., col=factor(s))) + fac_t + pal_s
p3 <- ggplot(fd, aes(x=sco_val, y=..scaled.., col=factor(t))) + fac_s + pal_t
p4 <- ggplot(fd, aes(x=sco_val, y=..scaled.., col=factor(s))) + fac_t + pal_s

ps <- list(
    (p1 / p2) & plot_annotation(tag_levels="a") & aes,
    (p3 / p4) & plot_annotation(tag_levels="a") & aes)

# saving
pdf(args[[2]], onefile=TRUE, width=15/2.54, height=12/2.54)
for (p in ps) print(p); dev.off()
