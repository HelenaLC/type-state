suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
df <- do.call(rbind, res)
# wrangling
.f <- \(df) df %>%
    group_by(dat, sel, sta) %>%
    summarise_at("sta_val", mean) %>%
    mutate(sta_val=case_when(
        sta_val < 0 ~ 0, 
        .default=sta_val))
df <- .f(df)

ps <- lapply(split(df, df$dat), \(fd) {
    ggplot(fd, aes(x = factor(sta), y = sta_val, color = factor(sel))) +
        scale_fill_brewer(palette="Dark2")+
        geom_bar(stat="identity", fill="white", position=position_dodge())+
        ggtitle(fd$dat[1]) + 
        xlab("Evaluation metrics") +
        ylab("Metric values") +
        labs(color="Feature selection")
})

pdf(args[[2]], onefile=TRUE, width=12, height=6)
for (p in ps) print(p); dev.off()