#!/usr/bin/env Rscript

# Self-contained clock-gene selection, analysis, data export, and plotting.
# The panel definition is intentionally local to this supplemental figure.

# Layout: code/ (this script), data/ (inputs; plotted numbers go to data/plotted/),
# figure/ (the images). Paths are relative to this figure's folder.
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
arg<-grep("^--file=",commandArgs(FALSE),value=TRUE);script<-normalizePath(gsub("~+~"," ",sub("^--file=","",arg),fixed=TRUE),mustWork=TRUE)
figdir<-dirname(script);root<-dirname(figdir)
dir.create(file.path(root,"data","plotted"),showWarnings=FALSE,recursive=TRUE);dir.create(file.path(root,"figure"),showWarnings=FALSE)
source(file.path(root,"code","analysis_universe.R"))
panel<-data.frame(
 gene_id=c("Sobic.007G047400","Sobic.004G216700","Sobic.003G040900","Sobic.001G411400","Sobic.002G275100","Sobic.005G044400","Sobic.003G191700","Sobic.001G340700","Sobic.009G137366","Sobic.001G352400","Sobic.003G166200"),
 label=c("LHY-like","TOC1-like","GIGANTEA","PRR73-like (PRR3/7 clade)","PRR95-like (PRR5/9 clade)","PRR59-like (PRR5/9 clade)","ELF3-like","ELF4-related","LUX/BOA-like","LNK1-like","LNK2-like"),stringsAsFactors=FALSE)
reference<-read.delim(gzfile(file.path(root,"data/generatable/external_diurnal_reference.tsv.gz")),check.names=FALSE)
metadata<-read.csv(file.path(root,"data/provided/sorghum_cold_timecourse_sample_metadata.csv"),check.names=FALSE)
current<-read.csv(gzfile(file.path(root,"data/externalsourcerequired/kallisto_tpms.csv.gz")),check.names=FALSE)
current_id<-sub("\\.[0-9]+$","",current$GeneID);reference<-reference[match(panel$gene_id,reference$normalized_gene_id),]
if(anyNA(reference$gene_id)|any(!panel$gene_id%in%current_id))stop("Clock genes failed to align.")
mat<-log2(as.matrix(current[match(panel$gene_id,current_id),-1,drop=FALSE])+PSEUDOCOUNT);storage.mode(mat)<-"numeric";rownames(mat)<-panel$gene_id;colnames(mat)<-names(current)[-1]
phasecols<-grep("^waveform_log2_tpm_zt",names(reference),value=TRUE);phases<-as.numeric(sub("waveform_log2_tpm_zt","",phasecols));wave<-as.matrix(reference[,phasecols]);rownames(wave)<-panel$gene_id
periodic<-function(values,p) periodic_waveform(values,phases,p)

# Confirm the locally selected LNK pair against the stated selection rule.
lnk_candidates<-c("Sobic.001G352400","Sobic.003G166200","Sobic.004G232100")
lnk_ref<-read.delim(gzfile(file.path(root,"data/generatable/external_diurnal_reference.tsv.gz")),check.names=FALSE)
lnk_ref<-lnk_ref[match(lnk_candidates,lnk_ref$normalized_gene_id),]
current_all<-as.matrix(current[,-1,drop=FALSE]);storage.mode(current_all)<-"numeric";rownames(current_all)<-current_id
control_mean<-rowMeans(current_all[lnk_candidates,metadata$Treatment=="CONT",drop=FALSE])
lnk<-data.frame(gene_id=lnk_candidates,current_control_mean_tpm=control_mean,
 external_mean_tpm=lnk_ref$mean_external_tpm,jtk_adjusted_p=lnk_ref$jtk_adjusted_p,
 jtk_period_h=lnk_ref$jtk_period_h,selected_for_figure=control_mean>=1&lnk_ref$jtk_adjusted_p<=.01&lnk_ref$jtk_period_h==24)
write.csv(lnk,file.path(root,"data","plotted","SFigure_clock_gene_lnk_selection.csv"),row.names=FALSE)
if(!identical(lnk$gene_id[lnk$selected_for_figure],tail(panel$gene_id,2)))stop("LNK selection changed.")

# Calibrate each external waveform to present chronological controls with one
# gene-wise slope/intercept; waveform shape and phase remain external.
control_meta<-metadata[metadata$Treatment=="CONT",];control_y<-mat[,metadata$Treatment=="CONT",drop=FALSE]
calibration<-data.frame(gene_id=panel$gene_id,slope=NA_real_,intercept=NA_real_)
for(i in seq_len(nrow(panel))){x<-vapply((control_meta$Time+1)%%24,function(p)periodic(wave[i,],p),numeric(1));fit<-lm(control_y[i,]~x);calibration$slope[i]<-coef(fit)[2];calibration$intercept[i]<-coef(fit)[1]}
write.csv(cbind(panel,reference[,c("mean_external_tpm","jtk_adjusted_p","jtk_period_h","jtk_lag_h","rhythmic_adjusted_p_le_0_01")],calibration[,-1]),file.path(root,"data","plotted","SFigure_clock_gene_panel.csv"),row.names=FALSE)

elapsed<-c(-1,3,7,11,15,19,23,27);rows<-list()
for(i in seq_len(nrow(panel)))for(onset in c("DAWN","DUSK"))for(e in elapsed){offset<-ifelse(onset=="DAWN",0,12);time<-if(e==-1)ifelse(onset=="DAWN",-1,11) else e+offset
 ctrl<-mat[i,metadata$Treatment=="CONT"&metadata$Time==time]
 cold<-if(e==-1)ctrl else mat[i,metadata$Treatment==onset&metadata$Time==time]
 zt<-(time+1)%%24;ext<-calibration$intercept[i]+calibration$slope[i]*periodic(wave[i,],zt)
 rows[[length(rows)+1]]<-data.frame(gene_id=panel$gene_id[i],label=panel$label[i],onset=onset,hours_after_onset=e,zt=zt,control_mean=mean(ctrl),cold_mean=mean(cold),external_prediction=ext)}
d<-do.call(rbind,rows);write.csv(d,file.path(root,"data","plotted","SFigure_clock_gene_trajectories.csv"),row.names=FALSE)

shade<-function(onset,ylim){for(a in seq(-1,27,by=.25)){zt<-(ifelse(onset=="DAWN",1,13)+a)%%24;rect(a,ylim[1],a+.25,ylim[2],border=NA,col=if(zt<14)adjustcolor("#F6D55C",.10)else adjustcolor("#4C78A8",.09))}}
draw<-function(device){device();on.exit(dev.off());layout(matrix(seq_len(22),nrow=11,byrow=TRUE));par(family="sans",mar=c(1.4,2.6,.9,.5),oma=c(2,0,0,0),mgp=c(1.7,.55,0),tcl=-.2,cex=.75)
 for(i in seq_len(nrow(panel)))for(onset in c("DAWN","DUSK")){z<-d[d$gene_id==panel$gene_id[i]&d$onset==onset,];yl<-range(c(z$control_mean,z$cold_mean,z$external_prediction),finite=TRUE);pad<-max(.15,diff(yl)*.08);yl<-yl+c(-pad,pad)
  plot(NA,xlim=c(-1,27),ylim=yl,xlab="",ylab="",xaxt="n",main="");shade(onset,yl);abline(h=pretty(yl),col=adjustcolor("grey70",.35),lwd=.4)
  lines(z$hours_after_onset,z$control_mean,col="#555555",lwd=1);points(z$hours_after_onset,z$control_mean,pch=1,cex=.65,col="#555555")
  col<-if(onset=="DAWN")"#D55E00" else "#0072B2";lines(z$hours_after_onset,z$cold_mean,col=col,lwd=1.3);points(z$hours_after_onset,z$cold_mean,pch=16,cex=.6,col=col)
  lines(z$hours_after_onset,z$external_prediction,col=col,lty=2,lwd=1)
  if(i==11)axis(1,at=elapsed,labels=elapsed) else axis(1,at=elapsed,labels=FALSE)
  if(onset=="DAWN")text(par("usr")[1]+.35,par("usr")[4]-.04*diff(par("usr")[3:4]),paste(panel$label[i],panel$gene_id[i],sep="   "),adj=c(0,1),cex=.58,font=2)
  if(i==1)mtext(onset,side=3,line=.1,font=2,cex=.8)
 }
 mtext("Hours after cold onset",side=1,outer=TRUE,line=.6,cex=.85)
}
draw(function()png(file.path(root,"figure","SFigure_clock_genes.png"),width=7.2,height=13,units="in",res=450,pointsize=9))
draw(function()svg(file.path(root,"figure","SFigure_clock_genes.svg"),width=7.2,height=13,pointsize=9))
