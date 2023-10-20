#args <- list(list.files("outs", "^sco-.*", full.names=TRUE), "plts/sco-ECDF.pdf")

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
        TRUE ~ sco_val, .default=sco_val))

df$sco_t <- paste(df$sco, df$t, sep="_")
df$sco_s <- paste(df$sco, df$s, sep="_")

# DE genes only (according to 'splatter')
de <- df[grep("GroupDE", names(df))]
fd <- df[!rowAlls(as.matrix(de) == 1), ]

# define 'true' type genes as
# having an average |logFC| > 1
# (using absolute values because we
# care about down-regulated markers)
is <- seq_len(ncol(de))
names(is) <- colnames(de)
fc <- vapply(is, \(i) 
    vapply(setdiff(is, i), \(j)
        log(de[, i]/de[, j], base=2), 
        numeric(nrow(de))) |> rowMeans(), 
    numeric(nrow(de)))
fc <- fd[rowMaxs(abs(fc)) > 1, ] 

gg <- list(
    labs(y="ECDF"),
    stat_ecdf(key_glyph="point"),
    scale_x_continuous(n.breaks=3),
    scale_y_continuous(n.breaks=3),
    facet_wrap(~ sco, scales="free_x", nrow=1),
    guides(color=guide_legend(override.aes=list(size=2))))

p1 <- ggplot(df, aes(sco_val, group=sco_t, col=factor(t))) 
p2 <- ggplot(fd, aes(sco_val, group=sco_t, col=factor(t))) 
p3 <- ggplot(fc, aes(sco_val, group=sco_t, col=factor(t))) 

p4 <- ggplot(df, aes(sco_val, group=sco_s, col=factor(s))) 
p5 <- ggplot(fd, aes(sco_val, group=sco_s, col=factor(s)))
p6 <- ggplot(fc, aes(sco_val, group=sco_s, col=factor(s)))

thm <- theme_linedraw(9) + theme(
    panel.grid=element_blank(),
    axis.title.x=element_blank(),
    legend.key.size=unit(0.5, "lines"),
    panel.spacing=unit(2, unit="mm"),
    strip.text=element_text(color="black", face="bold"),
    strip.background=element_rect(color=NA, fill="white"))

plt <- 
    wrap_elements(p1 / p2/ p3 + 
            plot_layout(guides="collect") & gg & thm & 
            scale_color_brewer(palette="Blues", "type\neffect")) / 
    wrap_elements(p4 / p5 / p6 + 
            plot_layout(guides="collect") & gg & thm & 
            scale_color_brewer(palette="Reds", "state\neffect")) + 
    plot_annotation(tag_levels="a") &
    theme(
        plot.margin=margin(0, unit="mm"),
        plot.tag=element_text(size=9, face="bold"))

ggsave(args[[2]], plt, units="cm", width=20, height=15)
