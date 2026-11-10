#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
script <- normalizePath(gsub("~+~"," ",sub("^--file=","",arg),fixed=TRUE),mustWork=TRUE)
# Layout: code/ (this script), data/ (inputs; plotted numbers go to data/plotted/),
# figure/ (the images). Paths are relative to this figure's folder.
analysis_dir <- dirname(dirname(script))   # this figure's folder
table_dir <- file.path(analysis_dir,"data","plotted"); figure_dir <- file.path(analysis_dir,"figure")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE); dir.create(figure_dir,showWarnings=FALSE)
source(file.path(analysis_dir,"code","analysis_universe.R"))

universe <- load_analysis_universe(analysis_dir,quiet=TRUE)
dawn <- rowMeans(response_matrix(universe,"DAWN",11))
# DUSK began 12 h later, so 11 h after DUSK onset was wall-clock T23.
dusk <- rowMeans(response_matrix(universe,"DUSK",23))

plot_data <- data.frame(
  gene_id=universe$genes,
  dawn_response_log2=dawn,
  dusk_response_log2=dusk,
  stringsAsFactors=FALSE
)
plot_data$direction <- ifelse(
  plot_data$dawn_response_log2 < 0 & plot_data$dusk_response_log2 > 0,
  "DAWN down / DUSK up",
  ifelse(plot_data$dawn_response_log2 > 0 & plot_data$dusk_response_log2 < 0,
         "DAWN up / DUSK down","Same direction")
)
direction_order <- c("Same direction","DAWN down / DUSK up","DAWN up / DUSK down")
plot_data$direction <- factor(plot_data$direction,levels=direction_order)
counts <- table(plot_data$direction)
# Category names stay as they are in the CSV; only the on-figure labels use Dawn/Dusk.
display_names <- c("Same direction","Dawn down / Dusk up","Dawn up / Dusk down")
legend_labels <- setNames(
  paste0(display_names," (",format(as.integer(counts[direction_order]),big.mark=","),")"),
  direction_order
)
r <- cor(plot_data$dawn_response_log2,plot_data$dusk_response_log2)

write.csv(plot_data,file.path(table_dir,"FigureS_11h_DAWN_DUSK_gene_responses.csv"),row.names=FALSE)
write.csv(data.frame(
  genes=nrow(plot_data),pearson_r=r,
  same_direction=as.integer(counts["Same direction"]),
  dawn_down_dusk_up=as.integer(counts["DAWN down / DUSK up"]),
  dawn_up_dusk_down=as.integer(counts["DAWN up / DUSK down"])
),file.path(table_dir,"FigureS_11h_DAWN_DUSK_summary.csv"),row.names=FALSE)

# Arial where installed, otherwise Liberation Sans (same letter widths).
font_family <- if (any(grepl("Arial", systemfonts::system_fonts()$family))) "Arial" else "Liberation Sans"
cols <- c("Same direction"="#B8B8B8",
          "DAWN down / DUSK up"="#008B72",
          "DAWN up / DUSK down"="#C51B7D")
p <- ggplot(plot_data,aes(dawn_response_log2,dusk_response_log2,color=direction))+
  geom_hline(yintercept=0,color="#555555",linewidth=.35)+
  geom_vline(xintercept=0,color="#555555",linewidth=.35)+
  geom_point(size=.65,alpha=.32)+
  geom_smooth(method="lm",formula=y~x,se=FALSE,color="black",linewidth=.75)+
  annotate("text",x=-7.5,y=-6.8,label=sprintf("Pearson r = %.3f\nn = %s",r,
    format(nrow(plot_data),big.mark=",")),hjust=0,size=3.8,family=font_family)+
  scale_color_manual(values=cols,breaks=direction_order,labels=legend_labels)+
  scale_x_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
  scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
  coord_fixed() +
  labs(
    x=expression("Dawn response at 11 h (log"[2]*" cold/control)"),
    y=expression("Dusk response at 11 h (log"[2]*" cold/control)"),
    color=NULL
  )+
  theme_classic(base_size=11,base_family=font_family)+
  theme(
    legend.position="top",
    legend.text=element_text(size=9),
    legend.key.width=unit(4,"mm"),
    axis.title=element_text(size=11),
    axis.text=element_text(size=10),
    plot.margin=margin(6,8,6,6)
  )+
  guides(color=guide_legend(nrow=1,override.aes=list(size=2.3,alpha=1)))

ggsave(file.path(figure_dir,"FigureS_11h_DAWN_DUSK_opposition.png"),p,
       width=6.5,height=6.2,dpi=300,bg="white")
# cairo_pdf embeds system fonts; the base pdf() device does not know Liberation Sans.
ggsave(file.path(figure_dir,"FigureS_11h_DAWN_DUSK_opposition.pdf"),p,
       width=6.5,height=6.2,bg="white",device=cairo_pdf)
svglite::svglite(file.path(figure_dir,"FigureS_11h_DAWN_DUSK_opposition.svg"),
                 width=6.5,height=6.2,pointsize=11)
print(p)
dev.off()

cat(sprintf("Wrote 11-h DAWN/DUSK opposition plot; r=%.9f; n=%d\n",r,nrow(plot_data)))
