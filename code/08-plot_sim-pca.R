args <- list(
    list.files("data/01-fil", full.names = TRUE)
)

lys <- lapply(args[[1]], readRDS)
