#!/usr/bin/env Rscript

# Supplemental diagnostic showing expression-dependent power to call rhythmicity
# and the rationale for the mean-external-TPM >= 3 analysis universe.

arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
script <- normalizePath(gsub("~+~", " ", sub("^--file=", "", arg), fixed=TRUE), mustWork=TRUE)
# Layout: code/ (this script), data/ (inputs; plotted numbers go to data/plotted/),
# figure/ (the images). Paths are relative to this figure's folder.
figdir <- dirname(script); root <- dirname(figdir)
dir.create(file.path(root, "data", "plotted"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root, "figure"), showWarnings=FALSE)
reference <- read.delim(
  gzfile(file.path(root,"data/generatable/external_diurnal_reference.tsv.gz")),
  check.names=FALSE, stringsAsFactors=FALSE
)
current <- read.csv(
  gzfile(file.path(root,"data/externalsourcerequired/kallisto_tpms.csv.gz")),
  check.names=FALSE, stringsAsFactors=FALSE
)
deg_results <- read.csv(
  gzfile(file.path(root,"data/generatable/deseq2_timepoint_contrasts.csv.gz")),
  check.names=FALSE, stringsAsFactors=FALSE
)
current_ids <- sub("\\.[0-9]+$","",current$GeneID)
deg_ids <- sub("\\.[0-9]+$","",deg_results$gene_id)
reference <- reference[
  reference$gene_id %in% intersect(current_ids, unique(deg_ids)),
]

# Equal-count bins reveal the detection curve without allowing the long TPM
# tail to determine bin widths.
reference <- reference[order(reference$mean_external_tpm,reference$gene_id),]
reference$expression_bin <- pmin(20L, ceiling(seq_len(nrow(reference))*20/nrow(reference)))
bins <- do.call(rbind,lapply(split(reference,reference$expression_bin),function(x) {
  n <- nrow(x); rhythmic <- sum(x$rhythmic_adjusted_p_le_0_01)
  # Wilson score interval (prop.test), which stays inside [0, 1] without clamping.
  interval <- prop.test(rhythmic, n)$conf.int
  data.frame(genes=n, median_mean_tpm=median(x$mean_external_tpm),
    rhythmic_fraction=rhythmic/n, lower_95=interval[1], upper_95=interval[2])
}))
thresholds <- c(.5,1,2,3,5,10)
threshold <- do.call(rbind,lapply(thresholds,function(v) {
  x <- reference[reference$mean_external_tpm>=v,]
  data.frame(minimum_mean_tpm=v,genes=nrow(x),rhythmic_genes=sum(x$rhythmic_adjusted_p_le_0_01),
    rhythmic_fraction=mean(x$rhythmic_adjusted_p_le_0_01))
}))
write.csv(bins,file.path(root,"data","plotted","SFigure_expression_rhythmicity_bins.csv"),row.names=FALSE)
write.csv(threshold,file.path(root,"data","plotted","SFigure_expression_rhythmicity_thresholds.csv"),row.names=FALSE)

draw <- function(device) {
  device(); on.exit(dev.off())
  par(mfrow=c(1,2),family="sans",mar=c(4.2,4.3,0.5,0.4),mgp=c(2.4,.7,0),tcl=-.25,cex=0.9)
  x <- bins$median_mean_tpm
  plot(x,100*bins$rhythmic_fraction,log="x",type="n",xlab="Mean expression in external series (TPM)",
    ylab="Genes identified as rhythmic (%)",ylim=c(0,100),xaxt="n")
  axis(1,at=c(.1,.3,1,3,10,30,100,300),labels=c("0.1","0.3","1","3","10","30","100","300"))
  polygon(c(x,rev(x)),100*c(bins$lower_95,rev(bins$upper_95)),border=NA,col=adjustcolor("#2C7FB8",.18))
  lines(x,100*bins$rhythmic_fraction,lwd=1.5,col="#2C7FB8"); points(x,100*bins$rhythmic_fraction,pch=16,cex=.75,col="#2C7FB8")
  abline(v=3,lty=2,col="#D55E00",lwd=1.2); text(3,5,"3 TPM cutoff",srt=90,pos=4,col="#D55E00",cex=.85)
  mtext("A",side=2,line=2.9,at=par("usr")[4],las=1,font=2,adj=1)

  plot(threshold$minimum_mean_tpm,100*threshold$rhythmic_fraction,type="b",pch=16,lwd=1.5,
    log="x",xlab="Minimum mean expression (TPM)",ylab="Retained genes identified as rhythmic (%)",
    ylim=c(75,90),xlim=c(.42,12),xaxt="n",col="#2C7FB8")
  axis(1,at=thresholds,labels=thresholds); abline(v=3,lty=2,col="#D55E00",lwd=1.2)
  text(threshold$minimum_mean_tpm,100*threshold$rhythmic_fraction,
    labels=format(threshold$genes,big.mark=","),pos=3,cex=.72,col="#444444")
  mtext("B",side=2,line=2.9,at=par("usr")[4],las=1,font=2,adj=1)
}
draw(function() png(file.path(root,"figure","SFigure_expression_rhythmicity.png"),width=7.2,height=3.3,units="in",res=600,pointsize=9))
draw(function() svg(file.path(root,"figure","SFigure_expression_rhythmicity.svg"),width=7.2,height=3.3,pointsize=9))
