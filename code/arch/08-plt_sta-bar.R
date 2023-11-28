#args <- list(list.files("outs", "^sta-", full.names=TRUE), "plts/sta-bar.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

idx <- grepl("-dat-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- do.call(rbind, res) %>%
    group_by(sel, sta, dat) %>%
    summarise_at("sta_val", mean, drop="none")

ps <- lapply(split(df, df$dat), \(fd)
    ggplot(fd, aes(x=sel, y=sta_val, col=sel)) +
        facet_wrap(~ sta, ncol=3, scales="free") + 
        geom_bar(stat="identity", fill="white") + 
        labs(title=fd$dat[1], x="selection method", y="evaluation values") +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)))

pdf(args[[2]], width=12, height=12, onefile=TRUE)
for (p in ps) print(p); dev.off()