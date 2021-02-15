rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)
if (!require("DESeq2")) BiocManager::install("DESeq2",update = T, ask=F)
if (!require("tximport")) BiocManager::install("tximport",update = T,ask=F)
if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap",update = T,ask=F)
if (!require("DEGreport")) {
  install.packages("lasso2",update=T)
  BiocManager::install("DEGreport",update = T,ask=F)}


if (!require(data.table)) install.packages("data.table")
if (!require("eulerr")) install.packages("eulerr")
if (!require("scales")) install.packages("scales")



library(DESeq2)
library(tximport)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(scales)

#Threshold for DEG
pthre<-0.01
FCthre<-1


path.out<-paste0("Output_alpha",pthre,"/")
if (dir.exists(path.out)) unlink(path.out, recursive=TRUE)
dir.create(path.out)

myoutput<-function(DESeqResult,file,tag,FCthre,pthre){
  res <- as.data.frame(DESeqResult)
  res<-cbind(data.frame(Gene=row.names(res)),res)
  write.table(res,
              file=file,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  
  res.sig<-res[(!is.na(res$padj)) & res$padj<pthre & abs(res$log2FoldChange)>FCthre,]
  file.sig<-paste0(dirname(file),"/Sig-",basename(file))
  write.table(res.sig,
              file=file.sig,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  
  
  res.sigmat<-data.table(Gene=rownames(res.sig),
                         Up=(res.sig$log2FoldChange >0),
                         Down=(res.sig$log2FoldChange <0)
  )
  names(res.sigmat)[2:3]<-paste(tag,names(res.sigmat)[2:3],sep=".")
  return(res.sigmat)
}

#Color Pallette for ggplot2
Lwd<-0.625
fontsize<-5
theme_pub<-theme_bw()+theme(
  aspect.ratio = 0.5,
  plot.margin = unit(c(5,5,5,5), "mm"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",size=1.5),
  text = element_text(face="bold",size = 3, colour = "black"),
  strip.background = element_blank(),
  strip.text=element_text(face="bold",size = fontsize, colour = "black"),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = fontsize, colour = "black"),
  axis.title.y.right = element_blank(),
  
  axis.ticks.y.left=element_line(),
  axis.ticks.y.right=element_line(),
  axis.ticks.length.y=unit(-5, "pt"),
  
  axis.text.x = element_text(size = fontsize, colour = "black",
                             hjust=0.5),
  axis.text.y = element_text(size = fontsize, colour = "black",
                             margin=unit(c(10,10,10,10), "pt")),
  axis.text.y.right = element_blank(),
  legend.position = "none",
  legend.title = element_text(size = fontsize),
  legend.text = element_text(size = fontsize),
  legend.text.align = 0
  #legend.position = "top"
)

pal<-c("black", # [1] WT/Mock
       "#FDD0A2", # [2] Mppif-KO #1
       "#FC8D59", # [3] gMpPIF/Mppif #1
       "#01665E", # [4] 35S:MpDELLA-Cit #1
       "#35978F", # [5] MpEF:MpDELLA #1
       "#80CDC1", # [6] MpEF:MpDELLA #5
       "#C7EAE5", # [7] MpEF:MpDELLA #6
       "#2171B5", # [8] gMpDELLA-GUS #1
       "#4292C6", # [9] gMpDELLA-GUS #2
       "#6BAED6", # [10] gMpDELLA-GUS #6
       "#3CA508", # [11] gMpDELLA-Cit
       "#009292"  # [12] DEX
       
)


#Get Mppif quant files
path.salmon<-"./Salmon_MpPIF"
file.quant.pifnTak<-list.files(path=path.salmon,pattern="\\.sf$",recursive = T,full.names = T)
file.quant.pif<-grep("pif",file.quant.pifnTak,value = T)
file.quant.Tak<-grep("Tak",file.quant.pifnTak,value = T)

#Generate tx2gene list
tx2gene.pif<-data.table(Transcript=fread(file.quant.pif[1])$Name)
tx2gene.pif[,Gene:=sub("(.*)\\..*","\\1",Transcript)]
tx2gene.Tak<-data.table(Transcript=fread(file.quant.Tak[1])$Name)
tx2gene.Tak[,Gene:=sub("(.*)\\..*","\\1",Transcript)]
tx2gene.common<-tx2gene.Tak[Gene %in% tx2gene.pif$Gene,]

# data table for making UpSet plot
UpSet.all<-data.table(Gene=unique(tx2gene.Tak$Gene))

#Mppif: Full interaction model
txi.pif.all<-tximport(file.quant.pif, type = "salmon", tx2gene = tx2gene.common)
txi.Tak.all<-tximport(file.quant.Tak, type = "salmon", tx2gene = tx2gene.common)
txi.all<-Map(cbind,txi.pif.all,txi.Tak.all)
txi.all$countsFromAbundance<-"no"
sampleTable.pif<-data.frame(Genotype = factor(rep(c("pif","Tak"),each=9),levels = c("Tak","pif")),
                        Time = factor(sub(".*FR([0-9]+)_.*","\\1",c(file.quant.pif,file.quant.Tak)))#,
                        # Sample = sub(".*\\/(.*)_zen.*","\\1",c(file.quant.pif,file.quant.Tak))
                        )

rownames(sampleTable.pif) <- colnames(txi.all$counts)
dds.all <- DESeqDataSetFromTximport(txi.all, sampleTable.pif, ~ Time+Genotype+Time:Genotype)
dds.all.keep<-rowSums(counts(dds.all)) >= 10
dds.all.kept <- dds.all[dds.all.keep,]
test.all <- DESeq(dds.all.kept)

#Visualization of sample distances
{
deseq2.rlg<-rlog(dds.all.kept,blind = F)
pcaData <- plotPCA(deseq2.rlg, intgroup=c("Genotype","Time"), returnData=TRUE)
pcaData$Batch<-paste0("Batch",sub(".*zen40([0-9]).*","\\1",c(file.quant.pif,file.quant.Tak)))
percentVar <- round(100 * attr(pcaData, "percentVar"))
p.pca<-ggplot(pcaData, aes(PC1, PC2, color=Time, shape=Genotype)) +
  geom_point(size=2*Lwd) +
  geom_text(aes(x=PC1+5,label=Batch),size=2,show.legend = F)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  # scale_color_manual(values = pal[c(1,2,3)]) +
  scale_shape_manual(values =c(1,16)) +
  coord_fixed() + theme_pub + 
  theme(aspect.ratio = 1,
        legend.position = "right",
        axis.title.x = element_text(size = fontsize, colour = "black"))
print(p.pca)
out.pdf<-paste0(path.out,"/Mppif_sample_distance_RLOG.pdf")
ggsave(file=out.pdf, plot=p.pca, width=10, height=10,units = "cm")
}

#Overall DEG by LRT:
{
  myClusterPlot<-function(degPatResult,ncol=5){
    cluster.df<-as.data.table(degPatResult[["normalized"]])
    cluster.df[,`:=`(cluster_n=length(unique(genes))),by=cluster]
    setorder(cluster.df,-cluster_n,cluster,genes)
    cluster.df[,`:=`(tag=paste("Group ",cluster,":",cluster_n," genes"))]
    cluster.df[,`:=`(tag=factor(tag,levels=unique(tag)),
                     Genotype=factor(Genotype,levels = c("Tak","pif")))]
    p.cluster<-ggplot(cluster.df,aes(x=Time, y=value, color = Genotype))  +
      geom_line(aes(group=paste0(genes,Genotype)),size=0.2*Lwd,alpha=0.3) +
      geom_smooth(aes(x = Time, y = value,group = Genotype, color = Genotype),
                  se = FALSE,alpha=0.5,size=0.7*Lwd,
                  method = "lm", formula = y~poly(x,2))+
      # geom_point(aes(group=Genotype),alpha = 0.3, size = Lwd,position = position_jitterdodge(dodge.width = 0.9)) +
      scale_color_manual(values = c("grey50",pal[2])) +
      geom_boxplot(data=cluster.df,aes(x=Time,y=value),alpha=0,size=0.5*Lwd)+
      scale_color_manual(values = pal[c(1,3)])+
      facet_wrap(~tag,ncol = ncol)
    
    p.cluster<-p.cluster+theme_pub+theme(aspect.ratio = 1,legend.position = "right")
    
    #Titles
    xtitle="Time"
    ytitle="Z-score if gene abundance"
    legendtitle = "Genotype"
    p.cluster<-p.cluster+labs(x=xtitle,y=ytitle,color=legendtitle)
    print(p.cluster)
    return(p.cluster)
  }
  
test.LRT.all<-DESeq(dds.all.kept,test="LRT",reduced = ~Time+Genotype)
res.LRT <- results(test.LRT.all,alpha = pthre)
sigmat.LRT<-myoutput(DESeqResult = res.LRT,
                     file = paste0(path.out,"Inoue_LRT.tsv"),
                     tag = "LRT",
                     pthre = pthre,
                     FCthre = 0
)

library(DEGreport)
rld_mat <- assay(deseq2.rlg)
cluster_rlog.pos <- rld_mat[rownames(rld_mat)%in%sigmat.LRT[LRT.Down==TRUE,Gene],]
clusters.pos <- degPatterns(cluster_rlog.pos, metadata = sampleTable.pif, time="Time", col="Genotype",
                        minc=0)
p.cluster.pos<-myClusterPlot(clusters.pos,ncol=5)
ggsave(filename = paste0(path.out,"Clusters-MpPIF_Pos.svg"),plot = p.cluster.pos,
       height = 20,width = 20,units = "cm")

cluster_rlog.neg <- rld_mat[rownames(rld_mat)%in%sigmat.LRT[LRT.Up==TRUE,Gene],]
clusters.neg <- degPatterns(cluster_rlog.neg, metadata = sampleTable.pif, time="Time", col="Genotype",
                            minc=0)
p.cluster.neg<-myClusterPlot(clusters.neg,ncol=5)
ggsave(filename = paste0(path.out,"Clusters-MpPIF_Neg.svg"),plot = p.cluster.neg,
       height = 20,width = 20,units = "cm")


}
#Helpful illustration on resultsNames & contrasts: https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
#Time-course comparison for each genotype
res.Tak1.FR1<-results(test.all, contrast = c("Time","1","0"),alpha = pthre)
sigmat.Tak1.FR1<-myoutput(DESeqResult = res.Tak1.FR1,
                             file = paste0(path.out,"Inoue_Tak1_FR1vs0.tsv"),
                             tag = "Tak1.FR1",
                             pthre = pthre,
                             FCthre = FCthre
                            )

res.Tak1.FR4<-results(test.all, contrast = c("Time","4","0"),alpha = pthre)
sigmat.Tak1.FR4<-myoutput(DESeqResult = res.Tak1.FR4,
                            file = paste0(path.out,"Inoue_Tak1_FR4vs0.tsv"),
                            tag = "Tak1.FR4",
                            pthre = pthre,
                            FCthre = FCthre
                          )
res.pif.FR1<-results(test.all, contrast = list( c("Time_1_vs_0","Time1.Genotypepif") ),alpha = pthre)
sigmat.pif.FR1<-myoutput(DESeqResult = res.pif.FR1,
                          file = paste0(path.out,"Inoue_pif_FR1vs0.tsv"),
                          tag = "pif.FR1",
                          pthre = pthre,
                          FCthre = FCthre
                          )
res.pif.FR4<-results(test.all, contrast = list( c("Time_4_vs_0","Time4.Genotypepif") ),alpha = pthre)
sigmat.pif.FR4<-myoutput(DESeqResult = res.pif.FR4,
                          file = paste0(path.out,"Inoue_pif_FR4vs0.tsv"),
                          tag = "pif.FR4",
                          pthre = pthre,
                          FCthre = FCthre
)


#Make UpSet plots to see the response of genes to FR
# Reference:https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
{
UpSet.FR<-Reduce(function(x,y) merge(x = x, y = y, by = "Gene",all.x=T), 
       list(UpSet.all,sigmat.Tak1.FR1,sigmat.Tak1.FR4,sigmat.pif.FR1,sigmat.pif.FR4))
UpSet.FR[is.na(UpSet.FR)] <- 0
UpSet.FR[,Sum:=rowSums(.SD, na.rm=T), .SDcols=names(UpSet.FR)[-1]]
UpSet.FR<-UpSet.FR[Sum>0,]

  cols.selected<-c("Tak1.FR1.Up","Tak1.FR4.Up",
                   "pif.FR1.Up","pif.FR4.Up",
                   "Tak1.FR1.Down","Tak1.FR4.Down",
                   "pif.FR1.Down","pif.FR4.Down")
  
  # m<-make_comb_mat(UpSet.pif[,..cols.selected],remove_complement_set = T)
  m<-make_comb_mat(UpSet.FR[,..cols.selected])
  UpSet(m,set_order = cols.selected)
  
  #Group by Up & Down
  t<-data.table(code=comb_name(m),code_n=as.numeric(comb_name(m)),cat=1,num=comb_size(m))
  t[,code:=factor(code,levels = code)]
  t[code_n>=10000,cat:=cat*2]
  t[(code_n %% 10000)>0,cat:=cat*3]
  setorder(t,cat,-num,-code_n)
  t[,od:=1:nrow(t)]
  setorder(t,code)
  m.order<-order(t$cat,-t$num,-t$code_n)
  
  #Make UpSet plot
  UpSet(m,set_order = cols.selected,
        comb_order = m.order)
  ss = set_size(m)
  cs = comb_size(m)
  ht = UpSet(m, 
             set_order = cols.selected,
             comb_order = m.order,
             width=unit(12, "cm"),
             height=unit(4, "cm"),
             pt_size=unit(3, "mm"),
             top_annotation = HeatmapAnnotation(
               "DEG Intersections" = anno_barplot(comb_size(m),
                                                  ylim = c(0, 420),
                                                  axis_param = list(
                                                    at = seq(0,400,100),
                                                    labels = seq(0,400,100),
                                                    labels_rot = 0),
                                                  border = TRUE,
                                                  gp = gpar(fill ="black"),
                                                  height = unit(4, "cm")
               ), 
               show_annotation_name = FALSE),
             left_annotation = rowAnnotation(
               " " = anno_text(ss, 
                               location = 1, 
                               just = "right",
                               width = max_text_width(ss) + unit(2, "mm")),
               "# of DEGs" = anno_barplot(-ss,
                                          baseline = 0,
                                          ylim = c(-1100,0),
                                          axis = FALSE,
                                          border = FALSE, 
                                          gp = gpar(fill = "black"),
                                          width = unit(2, "cm")
               ),
               set_name = anno_text(set_name(m), 
                                    location = 0.5, 
                                    just = "center",
                                    width = max_text_width(set_name(m)) + unit(4, "mm")),
               annotation_name_side = "top",
               annotation_name_rot = 0
             ), 
             right_annotation = NULL,
             show_row_names = F)
  pdf(paste0(path.out,"Mppif_Timecourse_Intersections.pdf"),height = 5,width =10 )
  ht = draw(ht)
  od = column_order(ht)
  
  # setorder(t2,od)
  decorate_annotation("DEG Intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = "centre", vjust=-0.5,
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
    grid.text("DEG Intersections", y = unit(1, "npc") + unit(2, "mm"), just = "bottom") 
  })
  dev.off()
  draw(ht)
}



#Mppif vs WT pair-wise comparison at each time point
res.FR0<-results(test.all, contrast = c("Genotype","pif","Tak"),alpha = pthre)
sigmat.FR0<-myoutput(DESeqResult = res.FR0,
                          file = paste0(path.out,"Inoue_FR0_pifvsTak.tsv"),
                          tag = "FR0",
                          pthre = pthre,
                          FCthre = FCthre
)
res.FR1<-results(test.all, contrast = list( c("Genotype_pif_vs_Tak","Time1.Genotypepif") ),alpha = pthre)
sigmat.FR1<-myoutput(DESeqResult = res.FR1,
                     file = paste0(path.out,"Inoue_FR1_pifvsTak.tsv"),
                     tag = "FR1",
                     pthre = pthre,
                     FCthre = FCthre
)
res.FR4<-results(test.all, contrast = list( c("Genotype_pif_vs_Tak","Time4.Genotypepif") ),alpha = pthre)
sigmat.FR4<-myoutput(DESeqResult = res.FR4,
                     file = paste0(path.out,"Inoue_FR4_pifvsTak.tsv"),
                     tag = "FR4",
                     pthre = pthre,
                     FCthre = FCthre
)

#Make UpSet plots to see the intersections at different time point
# Reference:https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
{
  UpSet.pif<-Reduce(function(x,y) merge(x = x, y = y, by = "Gene",all.x=T), 
                    list(UpSet.all,sigmat.FR0,sigmat.FR1,sigmat.FR4))
  UpSet.pif[is.na(UpSet.pif)] <- 0
  UpSet.pif[,Sum:=rowSums(.SD, na.rm=T), .SDcols=names(UpSet.pif)[-1]]
  UpSet.pif<-UpSet.pif[Sum>0,]
  UpSet.pif$FR.Res<-as.integer(UpSet.pif$Gene %in% UpSet.FR$Gene)
  
  
  cols.selected<-c("FR0.Up","FR1.Up","FR4.Up",
                   "FR0.Down","FR1.Down","FR4.Down")
  
  # m<-make_comb_mat(UpSet.pif[,..cols.selected],remove_complement_set = T)
  m<-make_comb_mat(UpSet.pif[,..cols.selected])
  UpSet(m,set_order = cols.selected)
  
  #Group by Up & Down
  t<-data.table(code=comb_name(m),code_n=as.numeric(comb_name(m)),cat=1,num=comb_size(m))
  t[,code:=factor(code,levels = code)]
  t[code_n>=1000,cat:=cat*2]
  t[(code_n %% 1000)>0,cat:=cat*3]
  t[cat==6,cat:=cat-(code_n%/%1000)%%10+code_n%%10]
  setorder(t,cat,-num,-code_n)
  t[,od:=1:nrow(t)]
  setorder(t,code)
  m.order<-order(t$cat,-t$num,-t$code_n)
  
  #Distinguish FR responsive genes
  cols<-c(cols.selected,"FR.Res")
  m1<-make_comb_mat(UpSet.pif[,..cols])
  t1<-data.table(code=comb_name(m1),num=comb_size(m1))
  t1[,`:=`(FR.Res=as.integer(substr(code,7,7)),
           code=substr(code,1,6))]
  
  t2<-data.table(dcast(t1, code~ FR.Res, value.var = "num",fill = 0))
  t2<-t2[code %in% comb_name(m),.(code,`1`,`0`)]
  t2<-merge(t2,t,by="code")
  t2[,code:=factor(code,levels = comb_name(m))]
  setorder(t2,code)
  
  
  #Make UpSet plot
  UpSet(m,set_order = cols.selected,
        comb_order = m.order)
  ss = set_size(m)
  cs = comb_size(m)
  ht = UpSet(m, 
             set_order = cols.selected,
             comb_order = m.order,
             width=unit(10, "cm"),
             height=unit(4, "cm"),
             pt_size=unit(3, "mm"),
             top_annotation = HeatmapAnnotation(
               "DEG Intersections" = anno_barplot(as.matrix(t2[,1:3],rownames = 1),
                                                  ylim = c(0, 350),
                                                  axis_param = list(
                                                    at = seq(0,300,100),
                                                    labels = seq(0,300,100),
                                                    labels_rot = 0),
                                                  border = TRUE,
                                                  gp = gpar(fill =c("black","grey50")),
                                                  height = unit(4, "cm")
               ), 
               show_annotation_name = FALSE),
             left_annotation = rowAnnotation(
               " " = anno_text(ss, 
                               location = 1, 
                               just = "right",
                               width = max_text_width(ss) + unit(2, "mm")),
               "# of DEGs" = anno_barplot(-ss,
                                          baseline = 0,
                                          ylim = c(-450,0),
                                          axis = FALSE,
                                          border = FALSE, 
                                          gp = gpar(fill = "black"),
                                          width = unit(2, "cm")
               ),
               set_name = anno_text(set_name(m), 
                                    location = 0.5, 
                                    just = "center",
                                    width = max_text_width(set_name(m)) + unit(4, "mm")),
               annotation_name_side = "top",
               annotation_name_rot = 0
             ), 
             right_annotation = NULL,
             show_row_names = F)
  pdf(paste0(path.out,"Mppif_Intersections.pdf"),height = 5,width =10 )
  ht = draw(ht)
  od = column_order(ht)
  
  decorate_annotation("DEG Intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = "centre", vjust=-0.5,
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
    grid.text("DEG Intersections", y = unit(1, "npc") + unit(2, "mm"), just = "bottom") 
  })
  dev.off()
  draw(ht)
}
{
  #Clusering of pif vs. WT FR4
  cluster_rlog.FR4.Up <- rld_mat[rownames(rld_mat)%in%sigmat.FR4[FR4.Up==TRUE,Gene],]
  clusters.FR4.Up <- degPatterns(cluster_rlog.FR4.Up, metadata = sampleTable.pif, time="Time", col="Genotype",
                                 minc=0)
  p.cluster.FR4.Up<-myClusterPlot(clusters.FR4.Up)
  ggsave(filename = paste0(path.out,"Clusters-Mppif_FR4_Up.svg"),plot = p.cluster.FR4.Up,
         height = 20,width = 20,units = "cm")
  cluster_rlog.FR4.Down <- rld_mat[rownames(rld_mat)%in%sigmat.FR4[FR4.Down==TRUE,Gene],]
  clusters.FR4.Down <- degPatterns(cluster_rlog.FR4.Down, metadata = sampleTable.pif, time="Time", col="Genotype",
                                   minc=0)
  p.cluster.FR4.Down<-myClusterPlot(clusters.FR4.Down)
  ggsave(filename = paste0(path.out,"Clusters-Mppif_FR4_Down.svg"),plot = p.cluster.FR4.Down,
         height =20,width = 20,units = "cm")
}


#Make a table for Mppif results
{
  file.list<-c("Inoue_FR0_pifvsTak.tsv",
               "Inoue_FR1_pifvsTak.tsv",
               "Inoue_FR4_pifvsTak.tsv",
               "Inoue_Tak1_FR1vs0.tsv",
               "Inoue_Tak1_FR4vs0.tsv",
               "Inoue_pif_FR1vs0.tsv",
               "Inoue_pif_FR4vs0.tsv"
  )
  file.list<-paste0(path.out,file.list)
  for (f in file.list){
    table<-fread(f)
    table.sele<-table[,.(Gene,log2FoldChange,padj)]
    names(table.sele)<-c("Gene",paste0(c("log2FC.","padj."),sub(".*Inoue_(.*).tsv","\\1",f)))
    if (f==file.list[1]) res.pif.all<-table.sele
    else res.pif.all<-merge(res.pif.all,table.sele)
  }
  write.table(res.pif.all,paste0(path.out,"Inoue_pif_all.tsv"),row.names = F,
              quote = F,
              sep = "\t")
}
#################### MpDELLAox #########################
file.IDs<-file.quant.Tak[1]
IDs.DELLA<-fread(file.IDs)

{
#Import data & DESeq analysis for DELLA
tx2gene.DELLA<-data.table(Transcript=IDs.DELLA$Name)
tx2gene.DELLA[,Gene:=sub("(.*)\\..*","\\1",Transcript)]

path.DELLA<-"./Salmon_MpDELLAox"
file.quant.DELLA<-list.files(path=path.DELLA,pattern="\\.sf$",recursive = T,full.names = T)

txi.DELLA <- tximport(file.quant.DELLA, type = "salmon", tx2gene = tx2gene.DELLA)
Genotags<-sub(".*\\/(.*)\\.quant\\.sf","\\1",file.quant.DELLA)

sampleTable.DELLA <- data.frame(Genotype = factor(sub("Mp_(.*)_.*","\\1",Genotags)))
rownames(sampleTable.DELLA) <- colnames(txi.DELLA$counts)

dds.DELLA <- DESeqDataSetFromTximport(txi.DELLA, sampleTable.DELLA, ~Genotype)
dds.DELLA.keep<-rowSums(counts(dds.DELLA)>= 1) >=3
dds.DELLA.kept<-dds.DELLA[dds.DELLA.keep,]

test.DELLA <- DESeq(dds.DELLA.kept)
res.DELLA<-results(test.DELLA, contrast = c("Genotype","OE","WT"),alpha = pthre)
sigmat.DELLA<-myoutput(DESeqResult = res.DELLA,
                       file = paste0(path.out,"MpDELLAox_OEvsWT.tsv"),
                       tag = "MpDELLAox",
                       pthre = pthre,
                       FCthre = FCthre
)


UpSet.both<-Reduce(function(x,y) merge(x = x, y = y, by = "Gene",all.x=T), 
                  list(UpSet.all,sigmat.LRT,UpSet.pif[,c("Gene","FR.Res","FR4.Up","FR4.Down")],sigmat.DELLA))
UpSet.both[is.na(UpSet.both)] <- 0
UpSet.both[,`:=`(pif.StrictUp=FR4.Up*FR.Res,
                 pif.StrictDown=FR4.Down*FR.Res
                 )]
}
{
#Venndiagram
library(eulerr)
myFisher<-function(A,B,v,total){
  # A<-"MpDELLAox.Up"
  # B<-"FR4.Up"
  x<-v$original.values[paste0(A,"&",B)]
  m<-sum(v$original.values[grep(A,names(v$original.values))])
  n<-total-m
  k<-sum(v$original.values[grep(B,names(v$original.values))])
  return(phyper(x-1,m,n,k,lower.tail = FALSE))
}

{
# Overlaps between MpDELLAox and Mppif vs WT (FR4)
select.cols<-c("MpDELLAox.Down","MpDELLAox.Up","FR4.Down","FR4.Up")
data.venn<-UpSet.both[,..select.cols]
v <- euler(data.frame(data.venn))
p.euler<-plot(v, fills = pal[c(7,5,3,2)],
        quantities = TRUE)
ggsave(paste0(path.out,"euler.svg"),plot = p.euler,width = 10,height = 10,units = "cm")

fishers<-data.table()
v.fisher<-euler(data.frame(UpSet.both[Gene %in% tx2gene.common$Gene,..select.cols]))
for (i in select.cols[1:2]){
  for (j in select.cols[3:4]){
    test<-myFisher(sub(" ",".",i),sub(" ",".",j),v.fisher,length(unique(tx2gene.common$Gene)))
    fishers<-rbind(fishers,data.table(set=names(test),
                                      value=v$original.values[match(names(test),names(v$original.values))],
                                      p=test))
  }
}
write.table(fishers,file = paste0(path.out,"euler-fishers.tsv"),row.names = F,quote = F,sep = "\t")
}
{
  # Overlaps between MpDELLAox and MpPIF targets
  select.cols<-c("MpDELLAox.Down","MpDELLAox.Up","LRT.Down","LRT.Up")
  data.venn.LRT<-UpSet.both[,..select.cols]
  v.LRT <- euler(data.frame(data.venn.LRT))
  p.LRT<-plot(v.LRT, fills = pal[c(7,5,3,2)],
          quantities = TRUE)
  ggsave(paste0(path.out,"euler.LRT.svg"),plot = p.LRT,width = 10,height = 10,units = "cm")
  
  fishers.LRT<-data.table()
  v.fisher.LRT<-euler(data.frame(UpSet.both[Gene %in% tx2gene.common$Gene,..select.cols]))
  
  for (i in select.cols[1:2]){
    for (j in select.cols[3:4]){
      test<-myFisher(sub(" ",".",i),sub(" ",".",j),v.fisher.LRT,length(unique(tx2gene.common$Gene)))
      fishers.LRT<-rbind(fishers.LRT,data.table(set=names(test),
                                        value=v$original.values[match(names(test),names(v$original.values))],
                                        p=test))
    }
  }
  write.table(fishers.LRT,file = paste0(path.out,"euler-fishers_LRT.tsv"),row.names = F,quote = F,sep = "\t")
}
}
