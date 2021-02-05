# R script to compare GO annotations in Marchantia
# Created 2020-01-22, by Rui Sun
# rm(list=ls())

if (!require(data.table,quietly =T)) install.packages("data.table")
if (!require(tidyr,quietly = T)) install.packages("tidyr")
if(!require(readxl)) install.packages("readxl")

library(data.table)

currDir<-getwd()
NothingButToDetectDirectory<-function(){}
setwd(getSrcDirectory(NothingButToDetectDirectory))

source("../00-Shared_Functions/GO_Operations.R")
source("../00-Shared_Functions/UpdateDatabase.R")


# update_go_default<-"n" # whether to update go.obo, "y" or "n" changes the default
# quiet_default<-T # whether to silent "myMessage"


# if (!exists("MyOntology")) {
#   load_go()
#   }

resultPATH<-"../../Blast2GO_Data/02-Results/CompareGO/"
b2gPATH<-paste0(resultPATH,"B2G/")
cleanPATH<-paste0(resultPATH,"B2G_clean/")
subPATH<-paste0(resultPATH,"Subtracted/")
statPATH<-paste0(resultPATH,"Stats/")

for (p in c(resultPATH,b2gPATH,cleanPATH,subPATH,statPATH)) {
  if (dir.exists(p)){
    unlink(sub("(.*)\\/$","\\1",p),recursive = T)
  }
  dir.create(p)
}

#Load JGI v3.1 genome and annotations
file.Genome<-"../../Blast2GO_Data/01-Parsed/Mp_genome_v3.1.tsv"
if (!file.exists(file.Genome)) parse_genome()
B2G.Genome<-fread(file.Genome)
write.table(B2G.Genome[!is.na(GO),.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"00-B2G_Genome.tsv"),
            sep = "\t",
            quote=F,row.names = F)

IDs<-unique(B2G.Genome[,.(Gene,Transcript,query)])

# Load blast2go annotations - Blast2GO:R
B2G.R<-fread("../../Blast2GO_Data/02-Results/Blast2GO/Blast2GO_clean.tsv")
B2G.R[,`:=`(Transcript=sub(".p", "", query,fixed = T),
            Gene=sub("\\..*", "", query))]
write.table(B2G.R[!GO=="",.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"01-B2G_R.tsv"),
            sep = "\t",
            quote=F,row.names = F)

# Load blast2go annotations - Interpro
B2G.Interpro<-fread("../../Blast2GO_Data/01-Parsed/Mp3.1_interpro.tsv")
# B2G.Interpro<-remove_redundancy(B2G.Interpro)
B2G.Interpro[,`:=`(Transcript=sub(".p", "", ID,fixed = T),
                   Gene=sub("\\..*", "", ID))]
write.table(B2G.Interpro[,.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"02-B2G_Interpro.tsv"),
            sep = "\t",
            quote=F,row.names = F)

# Load blast2go annotations - Blast2Go merged with Interpro
B2G.merged<-fread("../../Blast2GO_Data/02-Results/Blast2GO/Blast2GO+Interpro.tsv")
B2G.merged[,`:=`(Transcript=sub(".p", "", query,fixed = T),
                 Gene=sub("\\..*", "", query))]
write.table(B2G.merged[,.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"03-B2G_Merged.tsv"),
            sep = "\t",
            quote=F,row.names = F)

# Load blast2go annotations - Yamato2017, exported from .b2g
B2G.Yamato<-fread("../../Blast2GO_Data/01-Parsed/Marpol GO 20170815-exported.txt")
B2G.Yamato[,`:=`(Transcript=SeqName,
              Gene=sub("(.*)\\..*", "\\1", SeqName))]
B2G.Yamato.GO <- data.table(separate_rows(B2G.Yamato[,.(Gene,Transcript,`GO IDs`)],`GO IDs`,sep = "; "))
B2G.Yamato.GO[,GO:=sub(".*:(.*:)", "\\1", `GO IDs`)]
write.table(B2G.Yamato.GO[!GO=="",.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"04-B2G_Yamato.tsv"),
            sep = "\t",
            quote=F,row.names = F)


# Load blast2go annotations - Yamato2017, another copy
library(readxl)
B2G.Yamato1<-as.data.table(read_excel("../../Blast2GO_Data/01-Parsed/Marpol GO 20170815.xlsx"))
B2G.Yamato1[,`:=`(Transcript=`Sequence ID`,
                 Gene=sub("(.*)\\..*", "\\1", `Sequence ID`))]
B2G.Yamato1.GO <- data.table(separate_rows(B2G.Yamato1[,.(Gene,Transcript,GO)],GO,sep = "; "))
B2G.Yamato1.GO[,GO:=sub(".*:(.*:)", "\\1", GO)]
write.table(B2G.Yamato1.GO[!GO=="N/A",.(Gene,Transcript,GO)],
            file = paste0(b2gPATH,"05-B2G_Yamato1.tsv"),
            sep = "\t",
            quote=F,row.names = F)

# # Load blast2go annotations - PlantRegMap
# B2G.PRM<-fread("./Input/PlantRegMap_Mpo_GO_annotation")
# names(B2G.PRM)[1:2]<-c("Gene","GO")
# B2G.PRM.GO<-merge(B2G.PRM,IDs,by="Gene")
# write.table(B2G.PRM.GO[!GO=="",.(Gene,Transcript,GO)],
#             file = paste0(b2gPATH,"07-B2G_PRM.tsv"),
#             sep = "\t",
#             quote=F,row.names = F)

###############################################################################################

Sublist<-fread("./list-GA.csv",header = T)
setnames(Sublist,"tracking_id","Transcript")

Sublist.ont<-Sublist[rep(1:nrow(Sublist),3)]
Sublist.ont[,Ontology:=rep(c("biological_process", "cellular_component", "molecular_function"),each=nrow(Sublist))]

list<-list.files(b2gPATH)
for (file in list){
# for (file in list[1]){
  fileNAME<-paste0(b2gPATH,file)
  message(paste0("Processing ",file))
  Method=sub(".*_(.*)\\..*", "\\1", file)
  list.GO<-fread(fileNAME,header = T, sep = "\t")
  
  # Remove duplicated and redundant terms
  list.unique<-unique(fix_alternative(list.GO))
  list.clean<-remove_redundancy(list.unique,trSize = 1000,select.by = "Transcript")
  write.table(list.clean,
            file = paste0(cleanPATH,"Cleaned_",file),
            sep = "\t",
            quote=F,row.names = F)
  
  # Calculate redundancy for each transcript
  list.temp<-rbind(list.unique,list.clean)
  list.temp<-list.temp[,.(Red=.N==1),by=.(Gene,Transcript,GO)]
  list.summary<-list.temp[,.(Total=.N,Redundancy=sum(Red)),by=.(Gene,Transcript)]
  
  # Calculate specificity

  list.clean[,`:=`(NumAnc=num_anc(GO),
                   Ontology=get_prop(GO,"namespace"))]
  
  # Statistics
  stat<-data.frame(
    Method=Method,
    Duplication=1-nrow(list.unique)/nrow(list.GO),
    Redundancy=1-nrow(list.clean)/nrow(list.unique),
    Average_Redundancy=mean(list.summary[,Redundancy/Total]),
    Coverage=length(unique(list.unique$Gene))/length(unique(IDs$Gene)),
    Specificity_MF=mean(as.numeric(list.clean[Ontology=="molecular_function",NumAnc])),
    Specificity_BP=mean(as.numeric(list.clean[Ontology=="biological_process",NumAnc])),
    Specificity_CC=mean(as.numeric(list.clean[Ontology=="cellular_component",NumAnc]))
    
  )
  
  
  # Subsetting for genes of interest
  list.sub<-list.clean[Transcript %in% Sublist$Transcript]
  list.sub<-fix_obsolete(list.sub)
  list.sub[!is.na(GO),`:=`(
                 GO_Name=get_prop(GO,"name")
                 )]
  list.sub<-merge(Sublist.ont,list.sub,by=c("Ontology","Transcript"),all.x = T)
  list.sub[,Transcript:=factor(Transcript,levels = unique(Sublist$Transcript))]
  setkeyv(list.sub,c("Ontology","Transcript","GO"))
  write.table(list.sub[,.(Ontology,Transcript,name,GO,GO_Name)],
              file = paste0(subPATH,"Subtracted_",file),
              sep = "\t",
              quote=F,row.names = F)
  
  # Generate outputs for the subset
  sub<-list.sub[,.(Ontology,Transcript,name,GO,GO_Name)]
  sub[is.na(GO),GO_Name:=""]
  sub[,`:=`(GOs=paste0(GO," ",GO_Name),
             Method=Method)]
  sub<-sub[,.(GOs=paste0(GOs,collapse = "\n")),by=.(Method,Ontology,Transcript,name)]
  # write.table(sub,
  #             file = paste0(outPATH,"Combined_",file),
  #             sep = "\t",
  #             quote=T,row.names = F)

  if (file == list[1]) {
    stats<-stat
    specs<-list.clean[,Method:=Method]
    subs<-sub
  }else {
    stats<-rbind(stats,stat)
    specs<-rbind(specs,list.clean[,Method:=Method],fill=T)
    subs<-rbind(subs,sub)
    }
  
}

write.table(stats,
            file = paste0(statPATH,"CompareGO_stats.tsv"),
            sep = "\t",
            quote=F,row.names = F)
subs[,Method:=factor(Method,levels = sub(".*_(.*)\\..*", "\\1", list))]
setkeyv(subs,c("Method","Transcript","GOs"))
sub_wide <- spread(subs, Method, GOs)
write.table(sub_wide,
            file = paste0(subPATH,"CompareGO_sub.tsv"),
            sep = "\t",
            quote=T,row.names = F)


# setnames(specs,"Method","Trial")
library(ggplot2)
specs[,Method:=factor(Method,levels=unique(Method))]
vp<-ggplot(specs,aes(x=Method,y=NumAnc,fill=Ontology))+
  geom_violin()+facet_wrap(~Ontology,nrow = 1)+
  theme_bw()+
  theme(aspect.ratio = 0.8,
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1)
        )+
  scale_y_continuous(expand= c(0,0),limits=c(0,150),breaks = seq(0,150,30))+
  labs(x="Method",y="Specificity")
print(vp)
outFile.svg<-paste0(statPATH,"CompareGO_redundancies.svg")
ggsave(outFile.svg,plot = vp,width = 10,height = 5,units = "in")
# ggsave(outFile.pdf,plot = bp,width = 10,height = 10,units = "in")

setwd(currDir)
beep(3)


# #Extract certain GO terms
# input<-c("GO:0009684") #BP:indoleacetic acid biosynthetic process
# input<-c("GO:0009686") #BP:gibberellin biosynthetic process
# input<-c("GO:0009851") #auxin biosynthetic process


