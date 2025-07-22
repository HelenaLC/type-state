df <- readRDS("harmony_dat.rds")
library(ggplot2)
library(patchwork)

pal <- c(
    "#F4894B", 
    "#FAED6D",  
    "#D35FB7", 
    "#54C4DA"   
)

aes <- list(
    geom_bar(
        aes(fill=feature_selection, alpha=Harmony), key_glyph="point",
        stat="identity", position=position_dodge2(preserve="single")),
    guides(
        fill=guide_legend(override.aes=list(shape=21, stroke=NA, size=1)),
        alpha=guide_legend(override.aes=list(shape=21, stroke=NA, size=1))),
    scale_alpha_manual("Harmony\nintegration", values=c("No"=2/3, "Yes"=1)),
    scale_fill_manual(values=pal),
    labs(x=NULL, fill="selection"),
    theme_bw(4), theme(
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        plot.title=element_text(hjust=0.5)))

p1 <- ggplot(df, aes(
    reorder(feature_selection, LISI_g), LISI_g)) +
    ggtitle("mixing of cell\nstates/groups") + 
    labs(y="LISI_g (higher = better)") + aes

p2 <- ggplot(df, aes(
    reorder(feature_selection, -LISI_k), LISI_k)) +
    ggtitle("mixing of cell\ntypes/clusters") + 
    labs(y="LISI_k (lower = better)") + aes

gg <- wrap_plots(p2, p1, nrow=1) + 
    plot_layout(guides="collect") &
    scale_x_discrete(expand=expansion(0.2)) &
    scale_y_continuous(expand=expansion(0.03)) &
    theme(aspect.ratio=3/2,
        plot.background=element_blank(),
        legend.background=element_blank(),
        legend.spacing=unit(0, "pt"),
        legend.key.size=unit(0, "pt"),
        legend.key.spacing=unit(0, "pt"),
        legend.margin=margin(2, unit="pt"))
ggsave("~/projects/type-state/harmony.pdf", gg, units="cm", width=6, height=4)
