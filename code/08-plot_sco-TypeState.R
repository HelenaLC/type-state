suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
    library(stringr)
    library(wesanderson)
})

res <- lapply(args[[1]], \(x) { if (str_detect(x,"sim")) readRDS(x) })
res <- res[!vapply(res, is.null, logical(1))]
res <- do.call(rbind, res)

state <- res[res$sco == "state_edgeR", ] %>%
    rename('state_edgeR' = 'sco_val')
type <- res[res$sco == "type_Fstat", ] 

state$type_Fstat <- type$sco_val

## define type marker genes
gde <- state[grep("GroupDE", names(state))]
mg <- sapply(seq_len(ncol(gde)), \(i){
    not_i <- setdiff(seq_len(ncol(gde)), i)
    mgk <- sapply(not_i, \(j) {
        log(gde[,i]/gde[,j], base = 2)
    })
    rowMeans(mgk)
})

state$mg <- apply(abs(mg), 1, max)

## define state genes
cde <- state[grep("ConditionDE", names(state))]
cg <- sapply(seq_len(ncol(cde)), \(i){
    not_i <- setdiff(seq_len(ncol(cde)), i)
    cgk <- sapply(not_i, \(j) {
        log(cde[,i]/cde[,j], base = 2)
    })
    rowMeans(cgk)
})
state$cg <- apply(cg, 1, max)
pal <- wes_palette("Zissou1", 10000, type = "continuous")
gg <- ggplot(state, aes(type_Fstat, state_edgeR, col = cg )) +
    geom_point(shape = 16, alpha = 0.2, size = 0.7) + 
    facet_grid(t ~ s, labeller = \(.) label_both(.)) +
    #scale_colour_gradient2() +
    scale_colour_gradientn(colours = pal) + 
    scale_x_continuous() +
    scale_y_continuous() +
     theme_bw(9) + theme(
        #legend.position = "none",
        panel.grid = element_blank())

ggsave(args[[2]], gg, units = "cm", width = 15, height = 15)
