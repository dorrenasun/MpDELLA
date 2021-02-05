rm(list=ls(all.names = T))
if (!require(data.table,quietly =T)) install.packages("data.table")
if (!require(beepr,quietly =T)) install.packages("beepr")
if (!require(ontologyIndex,quietly =T)) install.packages("ontologyIndex")
if (!require(tidyr,quietly = T)) install.packages("tidyr")

library(data.table)
library(ontologyIndex)
library(tidyr)

source("../00-Shared_Functions/CommonUse.R")
source("../00-Shared_Functions/File_Locations.R")
source("../00-Shared_Functions/GO_Operations.R")
source("../00-Shared_Functions/UpdateDatabase.R")
##################################  Parameters  ##################################
ECw<-read.table("ECw.tsv",header = T,sep="\t",stringsAsFactors = F)
GOw<-5
args <- commandArgs(trailingOnly = TRUE)
# if() args <- c("--help")
if("--help" %in% args) 
{
  # cat()
  stop("Usage: Rscript Blast2GO-original.R --args parameterFile")
}else{
  while (length(args) < 1){
    args <- readline(prompt="Enter parameter file: ")
  }
  if (!file.exists(args)) stop("Parameter file does not exist!")
  print(args)
  Para<-fread(args)
  Blasttop<-Para[Name=="Blasttop",as.numeric(Value)]
  BestHit<-Para[Name=="BestHit",as.logical(Value)]
  Threshold<-Para[Name %in% c("cellular_component","biological_process","molecular_function"),]
  names(Threshold)<-c("Cat","Thres")
  Threshold[,Thres:=as.numeric(Thres)]
  cat(paste0("Blasttop: ",Blasttop,
             "\nBestHit: ",BestHit,
             "\nThreshold-CC: ",Threshold$Thres[1],
             "\nThreshold-BP: ",Threshold$Thres[2],
             "\nThreshold-MF: ",Threshold$Thres[3],
             "\n"
  ))
  print(Threshold)
}

update_databases<-"n" # whether to update go.obo, "y" or "n" changes the default
update_mapping<-"n" # whether to update mappings, "y" or "n" changes the default
update_score<-"n" # whether to update calculated scores, "y" or "n" changes the default

currDir<-getwd()
# dataPATH<-"../../Blast2GO_Data/"
outPATH<-path.blast2go
check_dir(outPATH)
BlastOpt<-paste0("_Blasttop",Blasttop,
                 ifelse(BestHit,"_BestHit",""))
outFile.blast<-paste0(outPATH,"Blast_best_hit",BlastOpt,".tsv")
outFile.mapping<-paste0(outPATH,"Blast2GO_mappings",BlastOpt,".tsv")
outFile.nomap<-paste0(outPATH,"Blast2GO_no_mapping",BlastOpt,".tsv")
outFile.score<-paste0(outPATH,"Blast2GO_scores",BlastOpt,".tsv")
outFile.final<-paste0(outPATH,"Blast2GO_clean",BlastOpt,
                      "_CC",Threshold$Thres[1],
                      "_BP",Threshold$Thres[2],
                      "_MF",Threshold$Thres[3],
                      ".tsv")
outFile.nonplant<-paste0(outPATH,"Blast2GO_non-plant",BlastOpt,
                         "_CC",Threshold$Thres[1],
                         "_BP",Threshold$Thres[2],
                         "_MF",Threshold$Thres[3],
                         ".tsv")
outFile.merged<-paste0(outPATH,"Blast2GO+Interpro_Blasttop",Blasttop,
                       "_CC",Threshold$Thres[1],
                       "_BP",Threshold$Thres[2],
                       "_MF",Threshold$Thres[3],
                       ".tsv")

# outPATH<-paste(outPATH,"/",sep = "")

file.blast<-file.Blast1.result.uniprot
file.uniprot.all<-file.uniprot.go
file.uniprot.filtered<-file.uniprot.go.selected

##################################  Updates  ##################################
update<-ifelse(exists("update_databases"),update_databases,NA)
while (! update %in% c("y","n")){
  update<-readline(prompt = "Update your uniprot database? (y/n)")
}
if (update=="y") {
  download_sprot()
  download_noiea()
  download_GO()
  download_taxonomy()
  update_taxonomy()
  update_uniprot() #depend on taxonomy and GO
  update_plantID() #depend on taxonomy and GO
  update_blastdb()
  
  download_genome()
  merge_genome_ver5_with_female()
  prepare_blast1()
}
if(!file.exists(file.uniprot.all)) update_uniprot()
if(!file.exists(file.uniprot.filtered)) update_plantID()


if (!file.exists(file.plantid)) update_plantID()

##################################  Main Program  ##################################
# GO.uniprot.all<-fread(file.uniprot.all,header = T)
GO.uniprot<-fread(file.uniprot.filtered,header = T)
ID.uniprot<-unique(GO.uniprot[,.(ID,Taxid)])
myMessage("Uniprot annotations loaded.")

update<-ifelse(exists("update_mapping"),update_mapping,NA)

if (!file.exists(outFile.mapping)) update<-"y"
while (! update %in% c("y","n")){
  update<-readline(prompt = "Update your mappings? (y/n)")
}

if (update=="y") {
  #Load GO annotations
  # GO.uniprot.plants<-fread(file.Dat.plants,header = T,sep = ",")
  
  #Load Blast results
  file.blast<-file.Blast1.result.uniprot
  if (!file.exists(file.blast)) stop("Blast file not found.")
  Blast.all<-read_blast(file.blast)
  # names(Blast.all)<-c("query", "ID", "pident", "ppos", "evalue", "bitscore")
  if (sum(!Blast.all$ID %in% GO.uniprot$ID)>0) stop("Blast results not match with mapping database.")
  
  Blast.all<-merge(Blast.all,ID.uniprot,by="ID",all.x = T)
  myMessage("Blast results loaded.")
  
  # Keep the best hit for each subject taxid
  Blast.1stHit<-Blast.all[!duplicated(Blast.all[,c("query","Taxid")]) & hit_num<=Blasttop,]
  write.table(Blast.1stHit,file=outFile.blast,quote = F,sep = "\t",row.names = F,col.names = F)
  myMessage("Best Blast hit from each species extracted.")

  # Mapping<-merge(Blast.1stHit,GO.uniprot,by=c("ID","Taxid"),allow.cartesian=T)
  if (BestHit) {
    Blast.selected<-Blast.1stHit
  }else {
    Blast.selected<-Blast.all[hit_num<=Blasttop,]
    }
  Mapping<-merge(Blast.selected,GO.uniprot,by=c("ID","Taxid"),allow.cartesian=T)
  Mapping<-fix_alternative(Mapping)
  write.table(Mapping,file=outFile.mapping,quote = F,sep = "\t",row.names = F)
  myMessage("Blast results mapped to uniprot annotation.")
  
  Blast.no_map<-Blast.selected[!query %in% Mapping$query,]
  write.table(Mapping,file=outFile.nomap,quote = F,sep = "\t",row.names = F)
  myMessage("Unmapped blast results reported.")
  
}else{
  Mapping<-fread(outFile.mapping)
  myMessage("Mappings loaded from file.")
}

update<-ifelse(exists("update_score"),update_score,NA)
if (!file.exists(outFile.score)) update<-"y"
while (! update %in% c("y","n")){
  update<-readline(prompt = "Update your scores? (y/n)")
}

if (update=="y") {
  myMessage("Calculating direct terms.")
  #Calculate DT
  Mapping[,DT:=ECw$Default[match(Evidence,ECw$EC)]*ppos]
  ECw.NAs<-Mapping[is.na(DT),unique(Evidence)]
  if (length(ECw.NAs)>0) stop("Evidence weight missing for: ",paste0(ECw.NAs,collapse=", "),". Please update'",file.Dat, "'.")
  
  #Add ancestral GO terms
  Mapping.all<-merge(Mapping,find_anc(Mapping$GO), by="GO",allow.cartesian=T)
  setkeyv(Mapping.all,c("query","Anc"))
  myMessage("Ancestral terms appended.")
  
  #Remove redundant entries for the same blast hit
  Mapping.unique<-Mapping.all[,.(Sim=max(ppos),
                                 Red=.N,
                                 DT=max(DT),
                                 EC=Evidence[which.max(DT)],
                                 Source=Source[which.max(DT)],
                                 Tag=paste0(unique(Tag),collapse="")
  ),by=.(query,ID,Anc)]
  setnames(Mapping.unique,"Anc","GO")
  myMessage("Redundant entries removed.")
  # Mapping.sub<-Mapping.unique[query=="Mapoly0097s0049.1.p"]
  
  #Find DT
  Score<-Mapping.unique[,.(NumHit=length(unique(ID)),
                           TopSim=ID[which.max(Sim)],
                           MaxSim=max(Sim),
                           TopDT=ID[which.max(DT)],
                           MaxDT=max(DT),
                           Tag=paste0(unique(Tag),collapse="")
  ),by=.(query,GO)]
  myMessage("Direct terms calculated.")
  # query.list<-query.list[1:220]
  
  # Calculate the numbers of subterms 
  num_child<-function(dt,trSize=100){
    # Check if GOs were loaded
    if (!exists("myOntology")) load_go()
    
    # dt<-Score
    # trSize<-5
    
    # Check for legal input
    if (sum(c("query","GO","Tag") %in% names(dt))<3) {
      stop("Check the input colomn names for 'query' and 'GO'.")
    }
    
    GO.list<-data.table(GO=unique(dt$GO))
    GO.list[,Des:=lapply(GO,function(x){paste0(get_descendants(myOntology,x),collapse = ",")})]
    GO.list<-data.table(separate_rows(GO.list,Des,sep = ","))
    GO.list<-GO.list[!Des==GO]
    GO.list.des<-GO.list[,.(Des=.(Des)),by=GO]
    
    TrueHits<-dt[Tag=="Ori",.(List=.(GO)),by=query]
    
    # for (i in 1:nrow(dt)){
    calculate_child<-function(q,go){
      set1<-unlist(GO.list.des[GO==go,Des])
      set2<-unlist(TrueHits[query==q,List])
      return(length(intersect(set1,set2)))
    }
    
    myMessage("Start calculating subnodes...")
    q<-0
    query.list<-unique(dt$query)
    dt$Subnode<-0
    setkey(dt,"query")
    while(q<length(query.list)){
      items<-which(dt$query %in% query.list[(q+1):min(q+trSize,length(query.list))])
      dt[items,Subnode:=mapply(calculate_child,query,GO)]
      q<-min(q+trSize,length(query.list))
      if (q %% 100==0 | q==length(query.list)) message(paste0(Sys.time(),"  Sub-nodes calculated (",q,"/",length(query.list),")."))
    }
    beep(3)
    myMessage("Subnodes calculated!")
    return(dt)
  }
  
  Score<-num_child(Score)
  setnames(Score,"Subnode","NumChild")
  # myMessage("Child entries calculated.")
  
  # Calculate final scores
  Score[,Score:=MaxDT+NumChild*GOw]
  Score[,Ont:=get_prop(GO,"namespace")]
  
  write.table(Score,file=outFile.score,quote = F,sep = "\t",row.names = F)
  myMessage("Final scores calculated.")
} else{
  Score<-fread(outFile.score)
  myMessage("Scores loaded from file.")
}

# Filter by category-specific thresholds
Score.filtered<-Score[Score>Threshold$Thres[match(Ont,Threshold$Cat)]]
Score.final<-remove_redundancy(Score.filtered)
setkeyv(Score.final,c("query","TopDT"))

write.table(Score.final,file=outFile.final,quote = F,sep = "\t",row.names = F)
message(paste0(Sys.time()," Results written to ",outFile.final))
stop()

#Find the annotations transferred from non-plant species.
{
source("../00-Shared_Functions/Tax_Operations.R")
if (!exists("myTaxon")) load_taxonomy()
plant_ids<-fread(file.plantid)
Score.final[,Taxid:=GO.uniprot$Taxid[match(TopDT,GO.uniprot$ID)]]
Score.nonplant<-Score.final[!Taxid %in% plant_ids$taxid]

Score.nonplant[,`:=`(Species=myTaxon$name[match(Taxid,myTaxon$taxid)],
                     GO_name=get_prop(GO,"name"))]

write.table(Score.nonplant,file=outFile.nonplant,quote = F,sep = "\t",row.names = F)
myMessage("Non-plant annotations reported.")
}

# Merge with Interpro results
file.interpro<-"../../Blast2GO_Data/01-Parsed/Mp3.1_interpro.tsv"
if (!file.exists(file.interpro)) stop("Interpro results not found!")
Interpro<-fread(file.interpro,header = T,sep = "\t")
setnames(Interpro,"ID","query")
Interpro[,query:=sub(".p","",query,fixed = T)]
Interpro<-fix_alternative(Interpro)

Merged<-unique(rbind(Score.final[,.(query,GO)],Interpro[,.(query,GO)]))
Merged.unique<-remove_redundancy(Merged,trSize = 1000)

write.table(Merged.unique,file=outFile.merged,quote = F,sep = "\t",row.names = F)
myMessage("Blast2GO Finished!")
setwd(currDir)

source("Blast2GO_statistics.R")


# source("../02-CompareGO/CompareGO.R")
beep(3)
