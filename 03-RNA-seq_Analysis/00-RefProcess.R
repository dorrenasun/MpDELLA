rm(list=ls())
if (!require(data.table,quietly =T)) install.packages("data.table")
if (!require(R.utils,quietly = T)) install.packages('R.utils')
library(R.utils)

library(data.table)

read_fasta<-function(file){
  
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  
  Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
              ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Text=paste(Full,collapse = "\n"),
                        Sequence=paste(Sequence,collapse = "")
  ),by=.(Num,ID)]
  return(Fasta.trans)
}

Dest<-"./RefGenome/"
if (!dir.exists(Dest)) dir.create(Dest)
  
file.v3<-paste0(Dest,"Mpolymorphav3.1.allTrs.fa.gz")
file.v3.genome<-paste0(Dest,"JGI_3.1.fasta.gz")
file.v5<-paste0(Dest,"MpTak1v5.1_r1.mrna.fasta")
file.v5.genome<-paste0(Dest,"MpTak1v5.1.fasta")

if (!file.exists(file.v3)) download.file(url = "https://marchantia.info/download/download/Mpolymorphav3.1.allTrs.fa.gz",destfile = file.v3)
if (!file.exists(file.v3.genome)) download.file(url = "https://marchantia.info/download/download/JGI_3.1.fasta.gz",destfile = file.v3.genome)

if (!file.exists(file.v5)) download.file(url = "https://marchantia.info/download/tak1v5.1/MpTak1v5.1_r1.mrna.fasta",destfile = file.v5)
if (!file.exists(file.v5.genome)) download.file(url = "https://marchantia.info/download/tak1v5.1/MpTak1v5.1.fasta",destfile = file.v5.genome)

fasta.v3<-read_fasta(file.v3)
fasta.v5<-read_fasta(file.v5)
fasta.v3.genome<-read_fasta(file.v3.genome)
fasta.v5.genome<-read_fasta(file.v5.genome)
#Save the V chromosome IDs from v5.1 to a file
write.table(grep("MpVg",fasta.v5$ID,value = T),file = paste0(Dest,"MpTak1_v5.1_V.txt"),col.names = F,row.names = F,quote = F)
#Remove V chromosom transcripts
fasta.v5.noV<-fasta.v5[-grep("MpVg",ID),]
#Save the X scaffold IDs from v3.1 to a file
list.X.scaff<-c("Mapoly0017s",
               "Mapoly0018s",
               "Mapoly0210s",
               "Mapoly0227s",
               "Mapoly0240s",
               "Mapoly0250s",
               "Mapoly0497s"
               )
list.X.Trs<-fasta.v3$ID[substr(fasta.v3$ID,1,11) %in% list.X.scaff]
fasta.v3.X<-fasta.v3[ID %in% list.X.Trs,]

fasta.Female<-rbind(fasta.v5.noV,fasta.v3.X)
write.table(fasta.Female$Text,file=paste0(Dest,"MpTak1_female.allTrs.fasta"),col.names=F,row.names=F,quote=F)

list.X.scaff.names<-c("scaffold_17",
                      "scaffold_18",
                      "scaffold_210",
                      "scaffold_227",
                      "scaffold_240",
                      "scaffold_250",
                      "scaffold_497")
fasta.Female.genome<-rbind(fasta.v5.genome[-grep("chrV",ID),],fasta.v3.genome[ID %in% list.X.scaff.names,])
write.table(fasta.Female.genome$Text,file=paste0(Dest,"MpTak1_female.genome.fasta"),col.names=F,row.names=F,quote=F)
fasta.Female.gentrome<-rbind(fasta.Female,fasta.Female.genome)
write.table(fasta.Female.gentrome$Text,file=paste0(Dest,"gentrome_female.fasta"),col.names=F,row.names=F,quote=F)
write.table(fasta.Female.genome$ID,file=paste0(Dest,"decoys_female.txt"),col.names=F,row.names=F,quote=F)

fasta.v5.gentrome<-rbind(fasta.v5,fasta.v5.genome)
write.table(fasta.v5.gentrome$Text,file=paste0(Dest,"gentrome_Tak1.fasta"),col.names=F,row.names=F,quote=F)
write.table(fasta.v5.genome$ID,file=paste0(Dest,"decoys_Tak1.txt"),col.names=F,row.names=F,quote=F)

fasta.all<-rbind(fasta.v5,fasta.v3.X)
fasta.all.genome<-rbind(fasta.v5.genome,fasta.v3.genome[ID %in% list.X.scaff.names,])
fasta.all.gentrome<-(rbind(fasta.all,fasta.all.genome))
write.table(fasta.all.gentrome$Text,file=paste0(Dest,"gentrome_all.fasta"),col.names=F,row.names=F,quote=F)
write.table(fasta.all.genome$ID,file=paste0(Dest,"decoys_all.txt"),col.names=F,row.names=F,quote=F)
write.table(fasta.all$ID,file=paste0(Dest,"Trascript_ids_all.txt"),col.names=F,row.names=F,quote=F)

