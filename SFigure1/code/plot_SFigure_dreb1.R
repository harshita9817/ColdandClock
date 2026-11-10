#!/usr/bin/env Rscript

# DREB1 (CBF) expression during chilling. DREB1 genes are nearly silent in
# unstressed leaves, so they fall outside the diurnal-adjusted analysis and are
# shown here from raw expression. Genes are the sorghum models whose closest
# Arabidopsis homolog is a DREB subfamily A-1 member and that pass the DESeq2
# testing filter.

# Layout: code/ (this script), data/ (inputs; plotted numbers go to data/plotted/),
# figure/ (the images). Paths are relative to this figure's folder.
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
arg<-grep("^--file=",commandArgs(FALSE),value=TRUE);script<-normalizePath(gsub("~+~"," ",sub("^--file=","",arg),fixed=TRUE),mustWork=TRUE)
figdir<-dirname(script);root<-dirname(figdir)
dir.create(file.path(root,"data","plotted"),showWarnings=FALSE,recursive=TRUE);dir.create(file.path(root,"figure"),showWarnings=FALSE)
PSEUDOCOUNT<-0.1
panel<-data.frame(
 gene_id=c("Sobic.002G269300","Sobic.002G269400","Sobic.002G269500","Sobic.002G269032","Sobic.007G181500","Sobic.010G016500"),
 label=rep("DREB1",6),stringsAsFactors=FALSE)
metadata<-read.csv(file.path(root,"data/provided/sorghum_cold_timecourse_sample_metadata.csv"),check.names=FALSE)
current<-read.csv(gzfile(file.path(root,"data/externalsourcerequired/kallisto_tpms.csv.gz")),check.names=FALSE)
if(!identical(names(current)[-1],metadata$`Sample name`))stop("TPM columns and metadata are not aligned.")
current_id<-sub("\\.[0-9]+$","",current$GeneID)
if(any(!panel$gene_id%in%current_id))stop("DREB1 genes missing from TPM matrix.")
tpm<-as.matrix(current[,-1]);storage.mode(tpm)<-"numeric"
mat<-log2(rowsum(tpm,current_id)[panel$gene_id,,drop=FALSE]+PSEUDOCOUNT)

elapsed<-c(-1,3,7,11,15,19,23,27);rows<-list()
for(i in seq_len(nrow(panel)))for(onset in c("DAWN","DUSK"))for(e in elapsed){offset<-ifelse(onset=="DAWN",0,12);time<-if(e==-1)ifelse(onset=="DAWN",-1,11) else e+offset
 ctrl<-mat[i,metadata$Treatment=="CONT"&metadata$Time==time]
 cold<-if(e==-1)ctrl else mat[i,metadata$Treatment==onset&metadata$Time==time]
 rows[[length(rows)+1]]<-data.frame(gene_id=panel$gene_id[i],label=panel$label[i],onset=onset,hours_after_onset=e,zt=(time+1)%%24,
  control_mean=mean(ctrl),control_se=sd(ctrl)/sqrt(length(ctrl)),cold_mean=mean(cold),cold_se=if(e==-1)sd(ctrl)/sqrt(length(ctrl)) else sd(cold)/sqrt(length(cold)))}
d<-do.call(rbind,rows);write.csv(d,file.path(root,"data","plotted","SFigure_dreb1_trajectories.csv"),row.names=FALSE)

shade<-function(onset,ylim){for(a in seq(-1,27,by=.25)){zt<-(ifelse(onset=="DAWN",1,13)+a)%%24;rect(a,ylim[1],a+.25,ylim[2],border=NA,col=if(zt<14)adjustcolor("#F6D55C",.10)else adjustcolor("#4C78A8",.09))}}
n<-nrow(panel)
draw<-function(device){device();on.exit(dev.off());layout(matrix(seq_len(2*n),nrow=n,byrow=TRUE));par(family="sans",mar=c(1.4,2.6,.9,.5),oma=c(2,1.2,0,0),mgp=c(1.7,.55,0),tcl=-.2,cex=.75)
 for(i in seq_len(n)){zz<-d[d$gene_id==panel$gene_id[i],];yl<-range(c(zz$control_mean-zz$control_se,zz$cold_mean+zz$cold_se,zz$cold_mean-zz$cold_se,zz$control_mean+zz$control_se),finite=TRUE);pad<-max(.15,diff(yl)*.08);yl<-yl+c(-pad,pad)
 for(onset in c("DAWN","DUSK")){z<-zz[zz$onset==onset,]
  plot(NA,xlim=c(-1,27),ylim=yl,xlab="",ylab="",xaxt="n",main="");shade(onset,yl);abline(h=pretty(yl),col=adjustcolor("grey70",.35),lwd=.4)
  col<-if(onset=="DAWN")"#D55E00" else "#0072B2"
  arrows(z$hours_after_onset,z$control_mean-z$control_se,z$hours_after_onset,z$control_mean+z$control_se,angle=90,code=3,length=.015,col="#555555",lwd=.6)
  arrows(z$hours_after_onset,z$cold_mean-z$cold_se,z$hours_after_onset,z$cold_mean+z$cold_se,angle=90,code=3,length=.015,col=col,lwd=.6)
  lines(z$hours_after_onset,z$control_mean,col="#555555",lwd=1);points(z$hours_after_onset,z$control_mean,pch=1,cex=.65,col="#555555")
  lines(z$hours_after_onset,z$cold_mean,col=col,lwd=1.3);points(z$hours_after_onset,z$cold_mean,pch=16,cex=.6,col=col)
  if(i==n)axis(1,at=elapsed,labels=elapsed) else axis(1,at=elapsed,labels=FALSE)
  if(onset=="DAWN")text(par("usr")[1]+.35,par("usr")[4]-.04*diff(par("usr")[3:4]),paste(panel$label[i],panel$gene_id[i],sep="   "),adj=c(0,1),cex=.62,font=2)
  if(i==1)mtext(onset,side=3,line=.1,font=2,cex=.8)
 }}
 mtext("Hours after cold onset",side=1,outer=TRUE,line=.6,cex=.85)
 mtext(expression(log[2]("TPM"+0.1)),side=2,outer=TRUE,line=0,cex=.85)
}
draw(function()png(file.path(root,"figure","SFigure_dreb1.png"),width=6.5,height=7.2,units="in",res=300,pointsize=9))
draw(function()svg(file.path(root,"figure","SFigure_dreb1.svg"),width=6.5,height=7.2,pointsize=9))
