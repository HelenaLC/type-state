# args <- list(
#     list.files("outs", "^sel-", full.names=TRUE),
#     "plts/sel-upset.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ComplexUpset)
})

idx <- grepl("-sim-", args[[1]])
res <- lapply(args[[1]][idx], readRDS)
df <- do.call(rbind, res) 

tsb <- c("t", "s", "b")
tsb <- lapply(tsb, \(.) sort(unique(df[[.]])))
tsb <- 100*do.call(expand.grid, tsb)
tsb <- apply(tsb, 1, \(.) sprintf("t%s,s%s,b%s", .[1], .[2], .[3]))
df$sim <- factor(df$sim, tsb)

.p <- \(df, 
    fill=NULL, sets=setdiff(names(df), fill), 
    main=NULL, n=20, size=0.8, lwd=0.2, font=6)
    upset(df, sets, name="", 
        set_sizes=FALSE, n_intersections=n, 
        height_ratio=1, min_degree=1, min_size=10, 
        themes=upset_default_themes(text=element_text(size=font)),
        base_annotations=list(
            "Intersection size"=
                intersection_size(counts=FALSE, mapping=aes(fill=.data[[fill]])) +
                ylab("") + ggtitle(main) +
                scale_fill_manual(values=c("lightgrey", "royalblue")) +
                theme(
                    legend.position="none", 
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_line(linewidth=0.2))),
        matrix=intersection_matrix(
            geom=geom_point(size=size), 
            segment=geom_segment(linewidth=lwd))) 

df$foo <- with(df, paste(sim, gene_id, sep=";"))
gs <- lapply(
    split(df, df$sel), 
    \(.) .$foo[.$sel_val])
tf <- UpSetR::fromList(gs)
tf$truth <- factor(tf$truth)
p1 <- .p(tf, fill="truth", n=80, size=2, lwd=0.5, font=9)
    
p2 <- lapply(split(df, df$sim), \(fd) {
    gs <- lapply(
        split(fd, fd$sel), 
        \(.) .$gene[.$sel_val])
    tf <- UpSetR::fromList(gs)
    tf$truth <- factor(tf$truth)
    .p(tf, fill="truth", main=sprintf("t%s,s%s", 100*fd$t[1], 100*fd$s[1]))
}) |> wrap_plots()

gg <- 
    wrap_elements(full=p1) / 
    wrap_elements(full=p2) + 
    plot_layout(heights=c(1, 3)) +
    plot_annotation(tag_levels="a") &
    theme(plot.tag=element_text(face="bold"))
ggsave(args[[2]], gg, width=30, height=35, units="cm")
