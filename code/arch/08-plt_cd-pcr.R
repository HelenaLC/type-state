#args <- list(list.files("data/01-fil", "cd\\.rd*", full.names=TRUE), "plts/cd-pcr.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggrastr)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
df <- bind_rows(res)

pcs <- grep("^PC", names(df), value=TRUE)
fig <- lapply(c("cluster_id", "group_id"), \(id) {
    gg <- lapply(split(df, df$sim), \(fd) {
        y <- as.matrix(fd[, pcs])
        fit <- summary(lm(y~fd[[id]]))
        R2 <- sapply(fit, \(.) .$adj.r.squared)
        data.frame(R2, PC=seq_along(pcs), t=fd$t[1], s=fd$s[1])
    }) |>
        do.call(what=rbind) |>
        #mutate(R2=case_when(R2 < 0 ~ abs(R2), .default=R2)) |>
        mutate(across(c("t", "s"), \(.) factor(., sort(unique(.))))) 
    
    ggplot(gg, aes(PC, 1, fill=R2)) +
        geom_tile(col="white", linewidth=0.2) +
        facet_grid(t ~ s, labeller=label_both) +
        scale_fill_gradientn(
            na.value="lightgrey",
            limits=c(0, 1), n.breaks=3,
            colors=rev(hcl.colors(9, "Blues"))) +
        scale_x_continuous(
            "principal component (PC)", 
            breaks=range(seq_along(pcs))) +
        coord_equal(length(pcs), expand=FALSE) + 
        ggtitle(sprintf("'%s' as predictor variable", id)) + 
        theme_bw(6) + theme(
            plot.background=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            panel.spacing=unit(0.25, "lines"),
            legend.key.height=unit(1, "lines"),
            legend.key.width=unit(0.5, "lines"),
            plot.title=element_text(hjust=0.5),
            strip.background=element_rect(fill="white"))
}) |> 
    wrap_plots(nrow=1) + 
    plot_layout(guide="collect") +
    plot_annotation(tag_levels="a") &
    theme(plot.tag=element_text(face="bold"))

ggsave("plts/cd-pcr.pdf", fig, width=15, height=8, units="cm")
