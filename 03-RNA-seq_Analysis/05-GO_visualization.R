rm(list = ls())

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (as.numeric(version$major)<4|(as.numeric(version$major)==4 & as.numeric(version$minor)<0.3)) {
  stop("Please update to R version 4.0.3 to make sure everything works properly.")}
BiocV<-as.numeric(paste(unlist(BiocManager::version()),collapse = "."))
if(BiocV<3.12) BiocManager::install(version = '3.12')
if(!require("GO.db")) BiocManager::install("GOSemSim",version = "3.12")

if (!require(GOSemSim)) BiocManager::install("GOSemSim",type="source",clean=T)

if (!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if (!require("AnnotationForge")) BiocManager::install("AnnotationForge",update=T,ask=F)
if (!require(data.table)) install.packages("data.table")

library(data.table)

source("Inputs/utilities.R") #A patch for the GOSimSem package
prepare_relation_df()
update<-F

library(AnnotationForge)
library(GO.db)
print(GO_dbInfo())

#Prepare annotation database for Marchantia
if (dir.exists("./org.Mpolymorpha.eg.db")) {
  unlink("./org.Mpolymorpha.eg.db",recursive = T)
}
if (file.exists("./org.Mpolymorpha.eg.sqlite"))   file.remove("./org.Mpolymorpha.eg.sqlite")

MpGO.ori<-fread("Inputs/Blast2GO_clean_RuiSun202011.tsv")
MpGO.ori[,gene:=sub("\\..*","",query)]

# Fix inconsistencies between the go.obo used during annotation and GO.db
ListSec<-data.table(GOid=keys(GO.db),Sec=Secondary(keys(GO.db)))
goids.NA<-data.table(go=unique(MpGO.ori[!GO %in% keys(GO.db),GO]))
goids.NA[go %in% keys(GOOBSOLETE), `:=`(Tag="obs",Ont=GOOBSOLETE[[go]]@Ontology),by=go]
goids.NA[!go %in% keys(GOOBSOLETE) & go %in% unlist(ListSec$Sec), `:=`(Tag="sec")]
if (sum(goids.NA$Tag=="sec",na.rm = T)>0) goids.NA[Tag=="sec",Ori:=ListSec[grep(pattern=go,Sec),GOid],by=go]
goids.NA[go=="GO:0042643",Ori:="GO:0005884"]

MpGO.ori[GO %in% goids.NA$go,GO:=goids.NA[go==GO,Ori],by=GO]

MpGO.sel<-unique(MpGO.ori[,.(gene,GO)])
MpGO.sel[,Evi:="IEA"]
names(MpGO.sel) <- c("GID","GO","EVIDENCE")

#Only GO annotations will be used, so there is no need to include genes with no GO annotations
MpSym<-MpGO.sel[,.(GID=unique(GID))]
MpSym[,`:=`(ENTREZID=GID,SYMBOL="-",GENENAME="-")]

MpChr<-MpGO.sel[,.(GID=unique(GID))]
MpChr[,CHROMOSOME:=ifelse(grepl("Mapoly",GID),"U",sub("Mp(.*)g.*","\\1",GID))]

if (update==T|!require(org.Mpolymorpha.eg.db))
{
  MpDb<-makeOrgPackage(gene_info=MpSym, chromosome=MpChr, go=MpGO.sel,
                       version="5.1.0",
                       maintainer="Rui Sun <dorrenasun@gmail.com>",
                       author="Rui Sun <dorrenasun@gmail.com>",
                       outputDir = ".",
                       tax_id="3197",
                       genus="Marchantia",
                       species="polymorpha",
                       goTable="go")
  install.packages(MpDb,repos = NULL,type="source")
  
}

  
if (!require("rrvgo")) BiocManager::install("rrvgo",update = F, ask=F)
if (!require("org.At.tair.db")) BiocManager::install("org.At.tair.db",update=T,ask=F)
if (!require("dendextend")) install.packages("dendextend")
if(!require(viridisLite))install.packages("viridisLite")
library(rrvgo)
library(dendextend)

library(viridisLite)

TopN<-25
orgdb<-"org.Mpolymorpha.eg.db"
method<-"Wang"
threshold<-0.77
Nclust<-12
outPath<-"./GO_vis/"
if (!dir.exists(outPath)) dir.create(outPath)

ReadTopGO<-function(inFile,Tag=""){
  # Tag<-"test"
  df<-fread(inFile,header = T)
  df[,`:=`(rank=1:.N,Tag=Tag)]
}
Overlap<-ReadTopGO("TopGO_a0.01/Overlap-FR4-Up_Fisher_BP.tsv",Tag="Up-Overlap")
DELLA<-ReadTopGO("TopGO_a0.01/MpDELLA-Up_Fisher_BP.tsv",Tag="MpDELLA-Up")
DELLA.1<-ReadTopGO("TopGO_a0.01/MpDELLA-Down_Fisher_BP.tsv",Tag="MpDELLA-Down")

PIF<-ReadTopGO("TopGO_a0.01/MpPIF_FR4-Up_Fisher_BP.tsv",Tag="Mppif-Up")
PIF.1<-ReadTopGO("TopGO_a0.01/MpPIF_FR4-Down_Fisher_BP.tsv",Tag="Mppif-Down")
GOTable<-do.call("rbind", list(DELLA,PIF,DELLA.1,PIF.1))
GOTable[,Tag:=factor(Tag,levels = unique(Tag))]
GOTable.red<-GOTable[rank<=TopN,]

{
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(svglite)) install.packages("svglite")
library(ggplot2)
library(ggrepel)
library(GOSemSim)
library(MASS) 
if (!require(Rtsne)) install.packages("Rtsne")
library(Rtsne)



MpGO<-godata(org.Mpolymorpha.eg.db,ont = "BP")
methods<-c("Resnik","Lin","Jiang","Rel","Wang")
  simMatrix <- calculateSimMatrix(unique(GOTable.red$GO),
                                  orgdb = orgdb,
                                  ont = "BP",
                                  semdata = MpGO,
                                  method=method
                                  )
  test<-as.dist(1-simMatrix)
{
    set.seed(3)
    tsne<-Rtsne(test,is_distance = T,perplexity = 5,theta=0.1)
    rownames(tsne$Y)<-rownames(simMatrix)
    
    GOpos<-data.table(GO.ID=rownames(tsne$Y),X=tsne$Y[,1],Y=tsne$Y[,2])
    
    hct<-as.dendrogram(hclust(dist(tsne$Y)))
    hct.cut<-dendextend:::cutree(hct,Nclust, order_clusters_as_data = F)
    cluster<-data.table(GO.ID=names(hct.cut),
                        Clust=hct.cut)
    hct<-color_branches(hct,Nclust,col=viridis(n=length(unique(cluster$Clust)),begin = 0,end=1),
                        groupLabels = T)

    plot(hct)
    
  Scores<-GOTable.red[,.(Term=unique(Term),Score=mean(-log10(classicFisher))),by=GO.ID]
  Scores<-merge(Scores,cluster,by="GO.ID")
  
  #Find topTerm
  listOff<-as.list(GOBPOFFSPRING)
  
  list.off<-listOff[Scores$GO.ID]  
  list.sum<-data.table(GO.ID=names(list.off),
                       Noff=sapply(list.off,length),
                       Ncover.own=sapply(1:length(list.off),function(i)sum(unlist(list.off[i]) %in% Scores[Clust==Clust[GO.ID==names(list.off)[i]],GO.ID])),
                       Ncover.All=sapply(list.off, function(x) sum(x %in% Scores$GO.ID))
  )
  Scores<-merge(list.sum,Scores,by="GO.ID")
  
  Scores[Ncover.own>1 & Ncover.own==Ncover.All,Tag:=1] #Specific terms
  Scores[Ncover.own>1 & Ncover.own<Ncover.All,Tag:=2] #Generic terms that cover other clusters
  Scores[Ncover.own<=1,Tag:=3] #Specific terms with no child
  Scores<-Scores[order(Clust,Tag,-Ncover.own,Noff,-Score),]
  
  Scores[,`:=`(topID=GO.ID[1],topTerm=Term[1]),by=Clust]
  Scores[,topTerm:=Term[match(topID,GO.ID)]]
  Scores<-Scores[order(Clust,-Score)]
  Scores[,`:=`(topTerm=factor(topTerm,levels = unique(topTerm)))]
  
  pdf(paste0(outPath,"Top",TopN,"_GO_clustering.pdf"),width = 10,height = 7)
  plot(hct)
  dev.off()
  
  GOpos.full<-merge(GOpos,Scores,by="GO.ID")
  
  df<-merge(GOTable.red,GOpos.full[,.(GO.ID,X,Y,topID,topTerm)],by="GO.ID")
  df[,Term:=factor(Term,levels = unique(Scores$Term))]

  #A plot showing all GO points  
  GOpos.cent<-GOpos.full[,.(X=mean(X),Y=mean(Y)),by=c("topTerm","Clust")]
  GOpos.full<-merge(GOpos,cluster,by="GO.ID")
  p0<-ggplot(GOpos.full,aes(x=X,y=Y))+
    geom_point(aes(color=Clust),alpha=0.6)+scale_color_viridis_c()+
    geom_text_repel(aes(label=topTerm),data=GOpos.cent)
  plot(p0)
  
}
# stop()



# Plotting


#Color Pallette
{
Lwd<-0.625
fontsize<-5
theme_pub<-theme_bw()+theme(
  aspect.ratio = 1,
  plot.margin = unit(c(5,5,5,5), "mm"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",size=1.5),
  text = element_text(face="bold",size = 3, colour = "black"),
  strip.background = element_blank(),
  strip.text=element_text(face="bold",size = fontsize, colour = "black"),
  # axis.title.x = element_blank(),
  axis.title = element_text(size = fontsize, colour = "black"),
  axis.title.y.right = element_blank(),
  
  axis.ticks.y.left=element_line(),
  axis.ticks.y.right=element_line(),
  axis.ticks.length=unit(-5, "pt"),
  
  axis.text.x = element_text(size = fontsize, colour = "black",
                             hjust=1, #angle = 45, 
                             margin=unit(c(10,10,10,10), "pt")),
  axis.text.y = element_text(size = fontsize, colour = "black",
                             margin=unit(c(10,10,10,10), "pt")),
  axis.text.y.right = element_blank(),
  # legend.position = "none",
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
}

#ReviGO-style plots
{
  p <- ggplot(df, aes(x=X, y=Y, color=topTerm)) + 
    geom_point(aes(size=-log10(classicFisher)), alpha=.5) +facet_wrap(~Tag,nrow = 2) +
    scale_color_manual(name="Cluster",values = viridis(length(unique(df$topTerm)), option = "D")) +
    scale_size_continuous(name="-log10(p)",range=c(1, 6)) +
    scale_x_continuous(n.breaks = 6,expand = expansion(mult=0.2)) +
    scale_y_continuous(n.breaks = 6,expand = expansion(mult=0.2)) +
    theme_pub+theme(legend.key.size = unit(2/length(unique(df$topTerm)),"cm"))  # p<-scatterPlot(simMatrix, reducedTerms)
  
  print(p)
  ggsave(filename = paste0(outPath,"Top",TopN,"_GO_Vis.svg"),plot = p,height = 10,width = 20,units = "cm")
  ggsave(filename = paste0(outPath,"Top",TopN,"_GO_Vis.png"),plot = p,height = 10,width = 20,units = "cm")
}

#Chart-style plots
{
  p1<-ggplot(df, aes(x=Tag, y=Term, size=Significant/Expected,color=-log10(classicFisher)))+
    geom_point()+
    scale_color_viridis_c(name="-log10(p)",alpha=0.8) +
    scale_size_continuous(range=c(2, 5)) +
    # scale_x_continuous(name="Semantic X",n.breaks = 6,expand = expansion(mult=0.2)) +
    # scale_y_continuous(name="Semantic Y",n.breaks = 6,expand = expansion(mult=0.2)) +
    theme_pub+theme(aspect.ratio = 1.5,axis.ticks.length.x=unit(5, "pt"))#,legend.position = "bottom")
  
  # print(p1)
  ggsave(filename = paste0(outPath,"Top",TopN,"_GO_Vis_chart.svg"),plot = p1,height = 20,width = 15,units = "cm")
  ggsave(filename = paste0(outPath,"Top",TopN,"_GO_Vis_chart.png"),plot = p1,height = 20,width = 15,units = "cm")
}
  
  
  cols<-c("Annotated","Significant","classicFisher","rank")
  GOTable.out<-dcast(GOTable[GO.ID %in% GO.ID[rank<=TopN] & classicFisher<0.05,], GO.ID ~ Tag,value.var =cols )
  new_col_order <- CJ(unique(GOTable.red$Tag), cols)[, paste(cols, V1, sep = "_")]
  
  GOTable.out<-merge(Scores[,.(GO.ID,Term,Clust,topTerm)],GOTable.out,by="GO.ID",sort=F)
  setcolorder(GOTable.out, c("Clust","topTerm","GO.ID","Term", new_col_order))
  setorderv(GOTable.out,c("Clust",grep("classicFisher",new_col_order,value = T)),na.last=T)
  write.table(GOTable.out,paste0(outPath,"Top_",TopN,"_enriched_GOs.tsv"),sep="\t",quote = F,row.names = F,col.names = T)
# p2<-ggplot(GOpos,aes(x=X,y=Y))+geom_point()
}
