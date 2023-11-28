#args <- list(list.files("outs", "^sco-.*", full.names=TRUE), "plts/sco-TypeState.pdf")

suppressPackageStartupMessages({
    library(dplyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
})

# loading
idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
res <- res[!vapply(res, is.null, logical(1))]
res <- do.call(rbind, res)

# wrangling
fd <- res |>
    filter(sco == "type_PVE")
df <- res |>
    filter(sco == "state_PVE") |>
    rename(state_PVE=sco_val) |>
    mutate(type_PVE=fd$sco_val)

de <- grep("^GroupDE", names(df))
ds <- grep("^ConditionDE", names(df))

df$de <- sapply(de, \(i) {
    j <- setdiff(de, i)
    lfc <- log2(df[[i]]/df[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

df$ds <- sapply(ds, \(i) {
    j <- setdiff(ds, i)
    lfc <- log2(df[[i]]/df[j])
    abs(rowMeans(lfc))
}) |> rowMeans()

# plotting
gg <- ggplot(df, aes(state_PVE, type_PVE, col=de)) +
    geom_point_rast(shape=16, alpha=0.2, size=0.7) + 
    facet_grid(t~s, labeller=\(.) label_both(.)) +
    scale_colour_gradientn(colors=hcl.colors(9, "Zissou1")) +
    theme_bw(9) + theme(panel.grid=element_blank())

# saving
ggsave(args[[2]], gg, units="cm", width=15, height=15)