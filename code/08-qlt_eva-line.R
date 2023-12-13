#args <- list(list.files("outs/dat", "^das-", full.names=TRUE), "plts/dat/das-upset.pdf")

suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)
df <- reshape2::melt(df, id = c("sel", "nGenes", "dat"))
names(df)[which(names(df)=="variable")] <- "sta"
names(df)[which(names(df)=="value")] <- "sta_val"


ps <- lapply(split(df, df$dat), \(fd) {
    ggplot(fd, aes(x = nGenes, y = sta_val, col = sel)) + 
        geom_line(stat = "identity", alpha=0.8) +
        geom_point(stat = "identity", alpha=0.8, size = 0.5) +
        facet_wrap(~sta, ncol=2, scales="free") +
        xlab("number of selected genes") + 
        theme_minimal() +
        scale_color_brewer(palette = "Set1") & theme(
            plot.margin=margin(),
            legend.position="bottom",
            legend.justification=c(0.5, 1),
            legend.box.spacing=unit(0, "pt"),
            panel.grid.minor=element_blank(), 
            plot.tag=element_text(size=9, face="bold"))
})



pdf(args[[2]], onefile=TRUE, width=10, height=8)
for (p in ps) print(p); dev.off()