if (!require(data.table,quietly =T)) install.packages("data.table")
if (!require(ggplot2,quietly =T)) install.packages("ggplot2")
if (!require(scales,quietly =T)) install.packages("scales")
if (!require(gridExtra,quietly =T)) install.packages("gridExtra")

library(data.table)
library(scales)
library(ggplot2)

# rm(list = ls())
currDir<-getwd()
NothingButToDetectDirectory<-function(){}
setwd(getSrcDirectory(NothingButToDetectDirectory))
if (!exists("read_fasta")) source("../00-Shared_Functions/UpdateDatabase.R")
source("../00-Shared_Functions/File_Locations.R")
source("../00-Shared_Functions/CommonUse.R")

##################################  Plot Statistics  ##################################
dataPATH<-"../../Blast2GO_Data/"
b2gPATH<-"../../Blast2GO_Data/02-Results/Blast2GO/"
outPATH<-"../../Blast2GO_Data/02-Results/Blast2GO_stats/"
check_dir(outPATH)

#Load total IDs from genome file
file.fasta<-file.fasta.query
if (!file.exists(file.fasta)) download_genome()
Total<-read_fasta(file.fasta)
Total[,Gene:=sub("\\..*", "", ID)]
setnames(Total,"ID","Transcript")

#Load Blast results
file.blast<-file.Blast1.result.uniprot
Blast<-read_blast(file.blast)
# names(Blast)<-c("query", "ID", "pident", "ppos", "evalue", "bitscore")
Blast[,Transcript:=sub(".p","",query,fixed = T)]
# setnames(Blast, "query_acc", "query")
Total[,Blast:=Transcript %in% Blast$Transcript]

#Load Mappings
outFile.mapping<-paste0(b2gPATH,"Blast2GO_mappings.tsv")
Mapping<-fread(outFile.mapping,header = T,sep = "\t")
Mapping[,Transcript:=sub(".p","",query,fixed = T)]
Total[,Mapping:=Transcript %in% Mapping$Transcript]

#Load Blast2GO annotations
outFile.blast2go<-paste0(b2gPATH,"Blast2GO_clean.tsv")
Blast2GO<-fread(outFile.blast2go,header = T,sep = "\t")
Blast2GO[,`:=`(Transcript=sub(".p", "", query,fixed = T),
            Gene=sub("\\..*", "", query))]

Total[,Blast2GO:=Transcript %in% Blast2GO$Transcript]

#Load Interpro results
file.Interpro<-"../../Blast2GO_Data/01-Parsed/Mp3.1_interpro.tsv"
Interpro<-fread(file.Interpro,header = T,sep = "\t")
Interpro[,Transcript:=sub(".p","",ID,fixed = T)]
Total[,Interpro:=Transcript %in% Interpro$Transcript]

# stop()
NoBlast<-Total[!(Blast|Interpro)==TRUE]#,.(Transcript,Sequence)]


Total.byGene<-Total[,.(Blast=sum(Blast)>0,
                       Mapping=sum(Mapping)>0,
                       Blast2GO=sum(Blast2GO)>0,
                       Interpro=sum(Interpro)>0
                       ),by=Gene]
# #Load Merged results
# file.merged<-paste0(outPATH,"Blast2GO+Interpro.tsv")
# Merged.unique<-fread(file.merged,header = T,sep = "\t")

#Summarize numbers of different categories
mySummary<-function(Total){
  Total.out<-Total[,.(Blast2GO=sum(Blast2GO),
                      `Mapping Only`=sum(Mapping & !Blast2GO),
                      `Blast Only`=sum(Blast & !Mapping),
                      Unknown=sum(!Blast),
                      Interpro=sum(Interpro),
                      `Blast2GO Only`=sum(Blast2GO & !Interpro),
                      Shared = sum(Blast2GO & Interpro),
                      `Interpro Only`=sum(Interpro & !Blast2GO),
                      `No Annotation`=sum(!Interpro & !Blast2GO),
                      Total=.N
  )]
  return(Total.out)
}
stats<-rbind(mySummary(Total),mySummary(Total.byGene))
stats$By<-c("Transcript","Gene")

stats.out<-data.table(t(stats[,-"By"]))
names(stats.out)<-stats$By
stats.out[,`:=`(Tag=names(stats)[1:length(names(stats))-1],
                `% Transcript`=percent(Transcript/stats[By=="Transcript",Total]),
                `% Gene`=percent(Gene/stats[By=="Gene",Total])
                )]
stats.out<-stats.out[,.(Tag,Transcript,`% Transcript`,Gene,`% Gene`)]

outFile.stats<-paste0(outPATH,"Blast2GO_stats.tsv")
write.table(stats.out,file=outFile.stats,quote = F,sep = "\t",row.names = F)

stats.out$Tag<-factor(stats.out$Tag,levels = rev(stats.out$Tag))


myPie<-function(stats,tag){
  stats$y<-stats[,..tag]
  stats[,`:=`(aes.x=ifelse(y/sum(y)<0.2,1.25,1),
              sum=cumsum(y),
              aes.y=y/2 + c(0, cumsum(y)[-length(y)]),
              label=percent(y/sum(y))
  )]
  label.color<-c("white",rep("black",nrow(stats[y!=0])-1))
  # if (length(label.color)>3) label.color[1]<-"white"
  
  bp<-ggplot(stats, aes(x=1, y=y, fill=Tag)) +
    geom_bar(width = 1, stat = "identity", color = "white")+
    geom_text(data=stats[y!=0],aes(x=aes.x,y = aes.y, label = label), color=label.color,size=4)+
    scale_y_continuous(expand= c(0,0),breaks = c(0,stats[y!=0]$sum))#+
  
  #theme for ggplot
  theme_pie<-theme_minimal()+theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(color="black", size=12, hjust=0.5, face="bold"),
    axis.text.y=element_blank(),
    axis.title = element_blank()
  )
  
  bp<-bp+theme_pie+coord_polar("y", start=0)
  return(bp)
}

test<-myPie(stats.out[1:4],"Gene")
test<-test+  labs(title="Blast2GO coverage\n(by transcript)") +scale_fill_brewer(palette = "Blues") 


# Piechart for Transcripts

bp.trans<-myPie(stats.out[1:4],"Transcript")+
  scale_fill_brewer(palette = "Blues") +
  labs(title="Blast2GO coverage\n(by transcript)")


# Piechart for Genes
bp.gene<-myPie(stats.out[1:4],"Gene")+
  scale_fill_brewer(palette = "Blues") +
  labs(title="Blast2GO coverage\n(by gene)")

# Piechart: comparing with Interpro, by transcript
bp.comp.trans<-myPie(stats.out[6:9],"Transcript")+
  scale_fill_brewer(palette = "Purples") +
  labs(title="Overlap with Interproscan\n(by transcript)")

# Piechart: comparing with Interpro, by gene
bp.comp.gene<-myPie(stats.out[6:9],"Gene")+
  scale_fill_brewer(palette = "Purples") +
  labs(title="Overlap with Interproscan\n(by gene)")

library(gridExtra)
bp<-grid.arrange(bp.trans, bp.gene,bp.comp.trans,bp.comp.gene, nrow = 2)

outFile.svg<-paste0(outPATH,"Blast2GO_summary.svg")
outFile.pdf<-paste0(outPATH,"Blast2GO_summary.pdf")
ggsave(outFile.svg,plot = bp,width = 10,height = 10,units = "in")
ggsave(outFile.pdf,plot = bp,width = 10,height = 10,units = "in")

myMessage("Blast2GO statistics calculated!")
setwd(currDir)

