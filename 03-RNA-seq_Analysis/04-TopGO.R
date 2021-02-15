rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

if (!require("topGO")) BiocManager::install("topGO",update = T, ask=F)
if (!require(data.table)) install.packages("data.table")
library(data.table)
library(topGO) # will also load GO.db

#Process GO annotations
MpGO.ori<-fread("Inputs/Blast2GO_clean_RuiSun202011.tsv")
MpGO.ori[,gene:=sub("\\..*","",query)]

# Fix inconsistencies between the go.obo used during annotation and GO.db
goids.NA<-data.table(go=unique(MpGO.ori[!GO %in% keys(GO.db),GO]))
goids.NA[go %in% keys(GOOBSOLETE), `:=`(Tag="obs",Ont=GOOBSOLETE[[go]]@Ontology),by=go]

ListSec<-data.table(GOid=keys(GO.db),Sec=Secondary(keys(GO.db)))
goids.NA[!go %in% keys(GOOBSOLETE) & go %in% unlist(ListSec$Sec), `:=`(Tag="sec")]
if (sum(goids.NA$Tag=="sec",na.rm = T)>0) goids.NA[Tag=="sec",Ori:=ListSec[grep(pattern=go,Sec),GOid],by=go]
goids.NA[go=="GO:0042643",Ori:="GO:0005884"]
MpGO.ori[GO %in% goids.NA$go,GO:=goids.NA[go==GO,Ori],by=GO]

MpGOlist<-MpGO.ori[,.(list=.(unique(GO))),by="gene"]
MpgeneID2GO<-MpGOlist$list
names(MpgeneID2GO)<-MpGOlist$gene

#Retrieve gene IDs
geneID.all <- fread("Inputs/Trascript_ids_all.txt",header = F)$V1
geneID.all<- unique(sub("\\..*","",geneID.all))
geneID.Tak1<-geneID.all[!grepl(pattern = "Mapoly",geneID.all)]
geneID.common<-geneID.Tak1[!grepl(pattern = "MpVg",geneID.Tak1)]

#Convert DESeq2 results to matrix
Res2Mat<-function(inFile,Tag="Data",alpha=0.01){
  DEGdata<-fread(inFile,header = T,sep = "\t")
  Sig.up<-DEGdata[!is.na(padj) & padj<alpha & log2FoldChange>1,Gene]
  Sig.down<-DEGdata[!is.na(padj) & padj<alpha & log2FoldChange<(-1),Gene]
  Mat<-DEGdata[,.(Gene,
                  Up=Gene %in% Sig.up,
                  Down=Gene %in% Sig.down
  )]
  names(Mat)<-c("Gene",paste0(Tag,c(".Up",".Down")))
  return(Mat)
}

pthre<-0.01
inDir<-paste0("Output_alpha",pthre,"/")

DELLA<-Res2Mat(paste0(inDir,"MpDELLAox_OEvsWT.tsv"),"DELLA",pthre)
PIF.LRT<-Res2Mat(paste0(inDir,"Inoue_LRT.tsv"),"LRT",pthre)
PIF.FR4<-Res2Mat(paste0(inDir,"Inoue_FR4_pifvsTak.tsv"),"FR4",pthre)

DEGMat<-data.table(Gene=geneID.all)
DEGMat<-Reduce(function(x,y) merge(x,y,by="Gene",all=T), list(DEGMat,DELLA,PIF.LRT,PIF.FR4))
DEGMat[is.na(DEGMat)]<-FALSE


#Classic Fisher enrichment analysis
myFisher<-function(ListGOI,geneID,Ont,topN=500,outFile){
  MpgeneList<-factor(as.integer(geneID %in% ListGOI))
  names(MpgeneList)<-geneID
  MpGOdata<-new("topGOdata",description="test",ontology=Ont,
                allGenes=MpgeneList,#geneSel = topDiffGenes,
                annot = annFUN.gene2GO,gene2GO=MpgeneID2GO,nodeSize = 5)
  resultFisher <- runTest(MpGOdata, algorithm = "classic", statistic = "fisher")
  resFisher<-as.data.table(GenTable(MpGOdata,classicFisher=resultFisher,orderBy="classicFisher",
                         topNodes=topN,numChar=100))
  write.table(resFisher,file=outFile,quote=F,row.names = F,sep="\t")
  outFile.simple<-paste0(dirname(outFile),"/Simplified_",basename(outFile))
  write.table(resFisher[,.(GO.ID,classicFisher)],file=outFile.simple,quote=F,sep="\t",row.names = F,col.names = F)
  return(resFisher)
}

{
outPath<-paste0("./TopGO_a",pthre,"/")
if(dir.exists(outPath)) unlink(sub("/$","",outPath),recursive = T)
dir.create(outPath)

res.DELLA<-myFisher(DEGMat[DELLA.Up|DELLA.Down,Gene],geneID.Tak1,
                       Ont = "BP",outFile = paste0(outPath,"MpDELLA_Fisher_BP.tsv"))
res.DELLA.Up<-myFisher(DEGMat[DELLA.Up==T,Gene],geneID.Tak1,
                          Ont = "BP",outFile = paste0(outPath,"MpDELLA-Up_Fisher_BP.tsv"))
res.DELLA.Down<-myFisher(DEGMat[DELLA.Down==T,Gene],geneID.Tak1,
                             Ont = "BP",outFile = paste0(outPath,"MpDELLA-Down_Fisher_BP.tsv"))
res.LRT<-myFisher(DEGMat[LRT.Up|LRT.Down,Gene],geneID.common,
                          Ont = "BP",outFile = paste0(outPath,"MpPIF_LRT_Fisher_BP.tsv"))
res.LRT.Up<-myFisher(DEGMat[LRT.Up==T,Gene],geneID.common,
                        Ont = "BP",outFile = paste0(outPath,"MpPIF_LRT-Up_Fisher_BP.tsv"))
res.LRT.Down<-myFisher(DEGMat[LRT.Down==T,Gene],geneID.common,
                           Ont = "BP",outFile = paste0(outPath,"MpPIF_LRT-Down_Fisher_BP.tsv"))
res.FR4<-myFisher(DEGMat[FR4.Up|FR4.Down,Gene],geneID.common,
                        Ont = "BP",outFile = paste0(outPath,"MpPIF_FR4_Fisher_BP.tsv"))
res.FR4.Up<-myFisher(DEGMat[FR4.Up==T,Gene],geneID.common,
                        Ont = "BP",outFile = paste0(outPath,"MpPIF_FR4-Up_Fisher_BP.tsv"))
res.FR4.Down<-myFisher(DEGMat[FR4.Down==T,Gene],geneID.common,
                           Ont = "BP",outFile = paste0(outPath,"MpPIF_FR4-Down_Fisher_BP.tsv"))
res.Ovrlap.FR4.Up<-myFisher(DEGMat[DELLA.Up & FR4.Up,Gene],geneID.common,
                           Ont = "BP",outFile = paste0(outPath,"Overlap-FR4-Up_Fisher_BP.tsv"))
res.Ovrlap.LRT.Up<-myFisher(DEGMat[DELLA.Up & LRT.Up,Gene],geneID.common,
                            Ont = "BP",outFile = paste0(outPath,"Overlap-LRT-Up_Fisher_BP.tsv"))
write.table(GO_dbInfo(),file=paste0(outPath,"GO_version_timestamp.txt"),row.names = F,sep = "\t",quote = F)
}
