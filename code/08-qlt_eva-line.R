#args <- list(list.files("outs/dat", "^eva-", full.names=TRUE), "plts/dat/eva-line.pdf")

suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(ggh4x)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)
df <- df %>% distinct(dat, num, sta, sel, sta_val, .keep_all = TRUE) %>%
  group_by(dat, num, sta, sel, nGenes) %>% 
  summarize(sta_val=mean(sta_val))


gg <- ggplot(df, aes(x = nGenes, y = sta_val, col = sel, group = sel)) + 
        geom_line(stat = "identity", alpha=0.8) +
        geom_point(stat = "identity", alpha=0.8, size = 0.5) +
        xlab("number of selected genes") + 
        theme_minimal() +
        ggh4x::facet_grid2(dat ~ sta, scales = "free", independent = "all") +
        scale_color_brewer(palette = "Paired") & theme(
            legend.position="bottom",
            legend.justification=c(0.5, 1),
            legend.box.spacing=unit(0, "pt"),
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA),
            plot.tag=element_text(size=9, face="bold"),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))



ggsave(args[[2]], gg, width=35, height=16, units="cm")