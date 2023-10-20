#args <- list(list.files("outs", "^sco-.*", full.names=TRUE), "plts/sco.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(matrixStats)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- bind_rows(res)
type_Fstat <- df %>%
    filter(sco == "type_Fstat") %>%
    pull("sco_val")
df <- df %>%
    filter(sco == "state_edgeR") %>%
    dplyr::rename(state_edgeR=sco_val) %>%
    mutate(type_Fstat)

# compute simulated 'type' & 'state' effects
# (average log2-fold change across clusters/groups)
.f <- \(df) {
    . <- seq_len(ncol(df))
    vapply(., \(i) 
        vapply(setdiff(., i), \(j)
            log2(df[, i]/df[, j]), 
            numeric(nrow(df))) |> 
            rowMeans(), 
        numeric(nrow(df)))
}
df$lfc_t <- rowMeans(abs(.f(df[grep("^GroupDE", names(df))])))
df$lfc_s <- rowMeans(abs(.f(df[grep("^ConditionDE", names(df))])))

ps <- lapply(paste0("lfc_", c("t", "s")), \(.)
    ggplot(df, aes(type_Fstat, state_edgeR, col=.data[[.]])) +
    facet_grid(t ~ s, labeller=\(.) label_both(.)) +
    geom_point(shape=16, alpha=0.2, size=0.7) + 
    scale_colour_gradientn(colours=hcl.colors(9, "Zissou 1")) + 
    theme_bw(9) + theme(panel.grid=element_blank()))

pdf(args[[2]], width=15/2.54, height=15/2.54, onefile=TRUE)
for (p in ps) print(p); dev.off()
