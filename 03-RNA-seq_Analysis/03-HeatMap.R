rm(list = ls())

if (!require(data.table)) install.packages("data.table")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(scales)) install.packages("scales")

library(data.table)
library(ggplot2)
library(scales)

{
  #Heatmap for Flavonoid genes
  file.list<-"Inputs/list-flavonoid-noPPO.csv"
  list<-fread(file.list)
  list[,`:=`(name=factor(name,levels = unique(name)),
             Gene=sub("\\..*","",ver5))]
  setorder(list,name,Gene)
  list[,text:=paste0(name," (",Gene,")")]
  list[,text:=factor(text,levels = rev(unique(text)))]
  
  pthre<-0.01
  FCthre<-1
  resdir<-paste0("./Output_alpha",pthre,"/")
  extractFC<-function(file, tag){
    res<-fread(paste0(resdir,file))
    res$Sample<-tag
    return(res[,.(Sample,Gene,log2FoldChange,padj)])
  }
  FC.DELLA<-extractFC("MpDELLAox_OEvsWT.tsv","MpDELLAox")
  FC.FR0<-extractFC("Inoue_FR0_pifvsTak.tsv","MppifvWT.FR0")
  FC.FR1<-extractFC("Inoue_FR1_pifvsTak.tsv","MppifvWT.FR1")
  FC.FR4<-extractFC("Inoue_FR4_pifvsTak.tsv","MppifvWT.FR4")
  FC.WT1<-extractFC("Inoue_Tak1_FR1vs0.tsv","WT.FR0")
  FC.WT4<-extractFC("Inoue_Tak1_FR4vs0.tsv","WT.FR4")
  FC.pif1<-extractFC("Inoue_pif_FR1vs0.tsv","Mppif.FR1")
  FC.pif4<-extractFC("Inoue_pif_FR4vs0.tsv","Mppif.FR4")
  
  data.all<-Reduce(function(x,y) rbind(x = x, y = y), 
                   list(FC.DELLA,FC.FR0,FC.FR1,FC.FR4))
  data.FR<-Reduce(function(x,y) rbind(x = x, y = y), 
                  list(FC.WT1,FC.WT4,FC.pif1,FC.pif4))
  data.all[,Sig:=(!is.na(padj))&(padj<pthre)&(abs(log2FoldChange)>FCthre)]
  data.FR[,Sig:=(!is.na(padj))&(padj<pthre)&(abs(log2FoldChange)>FCthre)]
  

  
  data.all[,FR.Res:=as.integer(data.all$Gene %in% data.FR$Gene[data.FR$Sig==T])]
  data.FR[,FR.Res:=as.integer(data.FR$Gene %in% data.all$Gene[data.all$Sig==T])]
  
  data.sub<-data.all[Gene %in% list$Gene]
  data.FR.sub<-data.FR[Gene %in% list$Gene]
  
  
  
  data.sub[,`:=`(text=list[match(data.sub$Gene,list$Gene),text],
                 name=list[match(data.sub$Gene,list$Gene),name],
                 Gene=factor(Gene,levels = (unique(list$Gene))),
                 Sample=factor(Sample,levels = c(rev(unique(Sample))))#,
                 # sig=(q_value<0.01)&(abs(`log2(fold_change)`)>1)
  )]
  data.FR.sub[,`:=`(text=list[match(data.FR.sub$Gene,list$Gene),text],
                 name=list[match(data.FR.sub$Gene,list$Gene),name],
                 Gene=factor(Gene,levels = (unique(list$Gene))),
                 Sample=factor(Sample,levels = c(rev(unique(Sample))))#,
                 # sig=(q_value<0.01)&(abs(`log2(fold_change)`)>1)
  )]

  Lwd<-0.625
  fontsize<-5
  
  #Theme for heatmap
  theme_hm<-theme_bw()+theme(
    plot.margin = unit(c(5,5,5,5), "mm"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    text = element_text(face="bold",size = 3, colour = "black"),
    strip.background = element_blank(),
    strip.text=element_text(face="bold",size = fontsize, colour = "black"),

    #Remove titles
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    #Remove x&y ticks
    axis.ticks = element_blank(),
    axis.ticks.length.x.top =unit(0, "pt"),
    axis.ticks.length.y =unit(0, "pt"),
    #Lengend position
    legend.position = "right",
    legend.key.size = unit(0.2, "cm"),
    #Format x and y texts
    axis.text.x.bottom =element_text(size = fontsize, angle = 45, hjust = 1,),
    axis.text.x.top =element_text(size = fontsize, angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    axis.text.y.left=element_text(face = "italic",size = fontsize, colour = "black"),
    axis.text.y.right=element_text(size = fontsize)
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
    
    
    myHM<-function(data){
      
      setorder(data,Gene,Sample)
      fir.labels<-unique(data$Gene)
      sec.labels<-data[match(fir.labels,Gene),name]
      sec.labels.dt<-data.table(Pos=1:length(sec.labels),Label=sec.labels)
      sec.labels.dt1<-sec.labels.dt[,.(Pos=mean(Pos),
                                       Start=min(Pos)-0.35,
                                       End=max(Pos)+0.35
      ),by=Label]
      
      data[,`:=`(y=as.integer(factor(Gene,levels = fir.labels)))]
      p<-ggplot(data, aes(x=y,y=Sample)) +
        geom_tile(aes(fill = log2FoldChange),color="black",width=0.9, height=0.9, size=0.3) +
        geom_text(data = data[Sig==T],aes(x=y,y=Sample),label="*",hjust=0.5,vjust=0.5)+
        geom_point(data = unique(data[,c("Gene","FR.Res","y")]),aes(x=y,y=0.2,color=factor(FR.Res,levels = c(0,1,2))),size=2*Lwd,show.legend=FALSE)+
        geom_segment(data=sec.labels.dt1,aes(x=Start,y = 4.6,xend=End,yend=4.6),color="black")+
        scale_fill_gradientn(colours = c(pal[3], "white", pal[5]),na.value="grey70",
                             breaks = c(-4,0,4), limits=c(-4.2,4.2),
                             oob=squish) +
        theme_hm+theme(aspect.ratio = (6)/length(unique(data$Gene)))
      
      p<-p+scale_y_discrete(expand = expansion(add = c(1, 0.7)))+
        scale_x_continuous(limits=c(0.4, length(fir.labels)+0.5),
                           expand = expansion(),
                           breaks = 1:length(fir.labels),labels = fir.labels,
                           sec.axis = sec_axis(~.,breaks = sec.labels.dt1$Pos,labels = sec.labels.dt1$Label))+
        scale_color_manual(breaks=c(0,1,2),values = c("grey70","black","orange"))
      
      return(p)
    }
    
    p.all<-myHM(data.sub)
    list.sig<-data.sub[,.(EverSig=sum(Sig)),by=.(Gene)]
    data.sub.sig<-data.sub[Gene %in% list.sig[(EverSig>0),Gene],]
    p.sig<-myHM(data.sub.sig)
    print(p.sig)
  
    p.FR.all<-myHM(data.FR.sub)
    # print(p.FR.all)
    list.FR.sig<-data.FR.sub[,.(EverSig=sum(Sig)+FR.Res),by=.(Gene)]
    
    p.FR.sig<-myHM(data.FR.sub[Gene %in% list.sig[(EverSig>0),Gene],])
    # print(p.FR.sig)

  
  # export svg files
  Fig.wd<-16
  Fig.ht<-4.5
  path.out<-"./HeatMap/"
  if(!dir.exists(path.out)) dir.create(path.out)
  tag<-sub(".*list-(.*)\\.csv","\\1",file.list)
  outFile.svg<-paste0(path.out,"Heatmap-",tag,"-all.svg")
  ggsave(file=outFile.svg, plot=p.all, width=Fig.wd, height=Fig.ht,units = "cm")
  ggsave(file=paste0(path.out,"Heatmap-",tag,"-all.png"), plot=p.all, width=Fig.wd, height=Fig.ht,units = "cm")
  outFile.svg<-paste0(path.out,"Heatmap-",tag,"-sig.svg")
  ggsave(file=outFile.svg, plot=p.sig, width=Fig.wd, height=Fig.ht,units = "cm")
  outFile.svg<-paste0(path.out,"Heatmap-",tag,"-FRsig.svg")
  ggsave(file=outFile.svg, plot=p.FR.sig, width=Fig.wd, height=Fig.ht,units = "cm")
}


