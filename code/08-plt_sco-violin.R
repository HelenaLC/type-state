#args <- list(list.files("outs", "^sco-.*", full.names=TRUE), "plts/sco.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

ex <- c("scmap")
df <- res |> bind_rows() |>
    mutate(sco_val=case_when(
        !(sco %in% ex) ~ log10(sco_val+1),
        .default=sco_val)) 

# DE genes only (according to 'Splatter')
de <- df[grep("GroupDE", names(df))]
fd <- df[!rowAlls(as.matrix(de) == 1), ]

gg <- list(
    geom_violin(trim=FALSE, fill="gray"),
    geom_boxplot(width=0.1),
    guides(color=guide_legend(override.aes=list(size=2))))

labs <- \(.) label_value(., multi_line=FALSE)
# facet_grid( ~ s, labeller=\(.) label_both(.), scales="free")
# facet_wrap(t ~ sco, labeller=labs, scales="free")
fac_s <- facet_wrap(s ~ sco, labeller=labs, scales="free") 
fac_t <- facet_wrap(t ~ sco, labeller=labs, scales="free") 

pal_s <- scale_color_brewer(palette="Reds", "state\neffect", limits=seq(0, 1, 0.2))
pal_t <- scale_color_brewer(palette="Blues", "type\neffect", limits=seq(0, 1, 0.2))

p1 <- ggplot(df, aes(x=factor(t), y=sco_val, col=factor(t))) + gg + fac_s + pal_t
p2 <- ggplot(fd, aes(x=factor(t), y=sco_val, col=factor(t))) + gg + fac_s + pal_t

p3 <- ggplot(df, aes(x=factor(s), y=sco_val, col=factor(s))) + gg + fac_t + pal_s
p4 <- ggplot(fd, aes(x=factor(s), y=sco_val, col=factor(s))) + gg + fac_t + pal_s

thm <- theme_bw(6) + theme(
    panel.grid=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_blank(),
    panel.spacing=unit(2, unit="mm"),
    plot.margin=margin(0, unit="mm"),
    legend.key.size=unit(0.5, "lines"),
    strip.text=element_text(color="black"),
    plot.tag=element_text(size=9, face="bold"),
    strip.background=element_rect(color=NA, fill="white"))

f1 <- (p1 / p3) & plot_annotation(tag_levels="a") & thm
f2 <- (p2 / p4) & plot_annotation(tag_levels="a") & thm

idx <- grepl("-dat-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- res |> bind_rows() |>
    filter(!is.na(sco_val)) %>%
    mutate(sco_val=case_when(
        !(sco %in% ex) ~ log10(sco_val+1),
        .default=sco_val)) 

ps <- lapply(split(df, df$dat), \(fd) 
    ggplot(fd, 
        aes(x=sco, y=sco_val)) + 
        gg + ggtitle(fd$dat[1]))

pdf(args[[2]], width=15/1.5, height=18/1.5, onefile=TRUE)
print(f1); print(f2); for (p in ps) print(p); dev.off()