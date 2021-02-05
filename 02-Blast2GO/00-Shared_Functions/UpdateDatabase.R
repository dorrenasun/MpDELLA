if (!require(zip,quietly = T)) install.packages('zip')
if (!require(R.utils,quietly = T)) install.packages('R.utils')
if (!require(data.table,quietly =T)) install.packages("data.table")
library(data.table)
if (!require(tidyr,quietly = T)) install.packages("tidyr")
library(tidyr)
if (!require(RCurl,quietly =T)) install.packages("RCurl")
if (!require(ggplot2,quietly =T)) install.packages("ggplot2")


NothingButToDetectDirectory<-function(){}
funDir<-paste0(getwd(),"/",getSrcDirectory(NothingButToDetectDirectory))
# funDir<-"../00-Shared_Functions"
print(funDir)

source("../00-Shared_Functions/CommonUse.R")
source("../00-Shared_Functions/File_Locations.R")


download_GO<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading gene ontologies.")
  check_dir(file.GO)
  download.file(url = "http://purl.obolibrary.org/obo/go.obo",destfile = file.GO)
  myMessage("Gene ontologies downloaded.")
  
  setwd(currDir)
}

download_genome<-function(){
  library(R.utils)
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading Marchantia genome.")
  check_dir(file.fasta.ver3)
  # file.zipped<-paste0(file.fasta.query,".gz")
  download.file(url = "http://marchantia.info/download/download/Mpolymorphav3.1.allTrs.pep_annot.fa.gz",
                # method = "wget",
                destfile = file.fasta.ver3)
  # if (file.exists(file.fasta.query)) file.remove(file.fasta.query)
  # gunzip(file.zipped)
  #
  check_dir(file.fasta.ver5)
  download.file(url = "https://marchantia.info/download/tak1v5.1/MpTak1v5.1_r1.protein.fasta",
                # method = "wget",
                destfile = file.fasta.ver5)

  myMessage("Marchantia genome downloaded.")
  
  setwd(currDir)
}

merge_genome_ver5_with_female<-function(){
  currDir<-getwd()
  setwd(funDir)
  myMessage("Merging Tak1 ver5.1 with female sequences.")
  
  if (!file.exists(file.fasta.ver3) | !file.exists(file.fasta.ver3)) download_genome()
  fasta.ver3<-read_fasta(file.fasta.ver3)
  fasta.ver5<-read_fasta(file.fasta.ver5)
  list.X.scaff<-c("Mapoly0017s",
                  "Mapoly0018s",
                  "Mapoly0210s",
                  "Mapoly0227s",
                  "Mapoly0240s",
                  "Mapoly0250s",
                  "Mapoly0497s"
  )
  list.X.Trs<-fasta.ver3$ID[substr(fasta.ver3$ID,1,11) %in% list.X.scaff]
  fasta.ver3.X<-fasta.ver3[ID %in% list.X.Trs,]
  fasta.ver5.wFe<-rbind(fasta.ver5, fasta.ver3.X)
  check_dir(file.fasta.ver5f)
  write.table(fasta.ver5.wFe$Text,file=file.fasta.ver5f,col.names=F,row.names=F,quote=F)
  
  myMessage("Tak1 ver5.1 merged with female sequences.")
  
  setwd(currDir)
}


parse_genome<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  if (!file.exists(file.fasta.query)) download_genome()
  Fasta<-fread(file.fasta.query,header = F,sep="\t",col.names = "Full")
  Total<-Fasta[grep(">",Full)]
  Total[,`:=`(Transcript=sub(" .*$","",sub(">","",Full)))]
  
  Total[,`:=`(Gene=sub("\\..*$","",Transcript),
              query=paste0(Transcript,".p"))]
  
  Total.GO<-data.table(separate_rows(Total[,.(Transcript,Full)],Full,sep="[,\\[\\]]"))
  Total.GO<-Total.GO[grep("GO:",Full,fixed = T)]
  setnames(Total.GO,"Full","GO")
  Total.final<-merge(Total[,Full:=NULL],Total.GO,all.x = T)
  
  check_dir(file.query.parsed)
  write.table(Total.final,file=file.query.parsed,quote = F,sep = "\t",row.names = F)
  myMessage("Marchantia genome Parsed.")
  setwd(currDir)
}



download_sprot<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading uniprot swissprot annotations.")
  check_dir(file.sprot.dat)
  download.file(url = "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
                destfile = file.sprot.dat)
  
  myMessage("Downloading swissprot fasta.")
  check_dir(file.sprot.fasta)
  download.file(url="https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                destfile = file.sprot.fasta)
  # download.file(url="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
  #               destfile = file.sprot.fasta)
  
  myMessage("Downloading swissprot date.")
  check_dir(file.sprot.reldate)
  download.file(url = "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt",
                destfile = file.sprot.reldate)
  
  myMessage("Swissprot data updated.")
  setwd(currDir)
}

download_noidea<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading all uniprot annotations without IEA.")
  check_dir(file.noiea.gaf)
  download.file(url="http://current.geneontology.org/annotations/goa_uniprot_all_noiea.gaf.gz",
                destfile = file.noiea.gaf)
  
  myMessage("Downloading release date for non-IEA.")
  check_dir(file.noiea.reldate)
  download.file(url="http://current.geneontology.org/summary.txt",
                destfile = file.noiea.reldate)
  
  myMessage("Uniprot non-IEA data updated.")
  setwd(currDir)
}

#Merge sprot and noiea GO annotations. Clean & fix alternative GO annotations
update_uniprot<-function(){
  {
  currDir<-getwd()
  setwd(funDir)
  if (!file.exists(file.sprot.fasta)) download_sprot()
  
  message("Release date for Swissprot:")
  date<-fread(file.sprot.reldate,sep = "",header = F)
  message(paste0(date$V1,"\n"))
  
  myMessage("Loading swissprot.")
  sprot.go<-read_genbank(file.sprot.dat,output = "GO")
  myMessage("Swissprot GO annotations loaded.")
  
  if (!file.exists(file.noiea.gaf)) download_noidea()
  
  message("Release date for non-IEA:")
  date<-fread(file.noiea.reldate,sep = "",header = F)
  message(paste0(date$V1,"\n"))
  
  myMessage("Loading uniprot - no IEA.")
  noiea.go<-read_gaf(file.noiea.gaf)
  myMessage("Non-IEA annnotations loaded.")
}
  # uniprot.go<-rbind(sprot.go,noiea.go)
  uniprot.go<-unique(rbind(sprot.go,noiea.go))
  myMessage("Uniprot annotations integrated.")
  
  if (!exists("fix_tax_IDs")) source("./Tax_Operations.R")
  uniprot.go[,Taxid:=fix_tax_IDs(Taxid)]
  uniprot.go<-unique(uniprot.go)

  if (!exists("fix_alternative")) source("./GO_Operations.R")
  uniprot.go<-fix_alternative(uniprot.go)
  
  uniprot.clean<-remove_redundancy(uniprot.go,trSize=5000,select.by = c("ID","Evidence"))
  myMessage("GO annotations cleaned.")
  
  check_dir(file.uniprot.go)
  write.table(uniprot.clean,file=file.uniprot.go,quote = F,sep = "\t",row.names = F)
  myMessage("Uniprot database updated.")
  
  setwd(currDir)
}


update_plantID<-function(){
  {
  currDir<-getwd()
  setwd(funDir)
  if (!file.exists(file.uniprot.go)) update_uniprot()
  uniprot.all<-fread(file.uniprot.go)
  
  if (!exists("load_taxonomy")) source("Tax_Operations.R")
  if (!exists("myTaxon")) load_taxonomy()
  uniprot.all[,Taxid:=fix_tax_IDs(Taxid,sum_at = "species")]
  
  list.unknown<-unique(uniprot.all[!Taxid %in% myTaxon$taxid,.(Taxid,ID)])
  
  
  if (file.exists(file.unknown)) {
    list.unknown.ex<-fread(file.unknown)
    list.unknown.new<-list.unknown[!ID %in% list.unknown.ex$ID]
    list.unknown<-rbind(list.unknown.ex,list.unknown.new[,is_plant:=""])
  } else {
    list.unknown[,is_plant:=""]
  }
  check_dir(file.unknown)
  write.table(list.unknown,file=file.unknown,quote = F,sep = "\t",row.names = F)
  # message("Please check Blast2GO_Data/01-Parsed/Taxonmy/Manual_taxids.tsv for unkown taxids.")
  
  list.manual<-fread(file.unknown)
  while (sum(is.na(list.manual$is_plant))>0) {
    message("Unknown taxids needs manual lookup in Blast2GO_Data/01-Parsed/Taxonmy/Manual_taxids.tsv.")
    invisible(readline(prompt="Press any key to continue."))
    list.manual<-fread(file.unknown)
  }
  
  plant_ids<-get_tax_descendants(33090)  #green plants
  
  plant_ids.all<-data.table(Taxid=union(plant_ids$taxid,list.manual[is_plant==TRUE,Taxid]))
  }
  check_dir(file.plantid)
  write.table(plant_ids.all,file.plantid,quote=F,row.names = F)
  
  if (!exists("get_prop")) source("GO_Operations.R")
  uniprot.all[,Cat:=get_prop(GO,"namespace")]
  list.GO.plant.BP<-find_anc(uniprot.all[Cat=="biological_process" & Taxid %in% union(plant_ids$taxid,list.manual[is_plant==TRUE,Taxid]),GO])
  uniprot.all[GO %in% list.GO.plant.BP$Anc,Cat:="is_plant"]
  list.GO.nonplant<-unique(uniprot.all[Cat=="biological_process",.(Taxid,GO)])
  list.GO.nonplant[,TaxName:=paste0(myTaxon[match(Taxid,taxid),name]," (",Taxid,")")]
  list.GO.nonplant<-list.GO.nonplant[,.(Taxons=paste(TaxName,collapse = ", ")),by="GO"]
  list.GO.nonplant[,GO_Name:=get_prop(GO,"name")]
  check_dir(file.nonplant.BP)
  write.table(list.GO.nonplant, file=file.nonplant.BP,quote=F,sep = "\t",row.names = F)
  
  
  uniprot.filtered<-uniprot.all[!Cat=="biological_process",.(ID,Taxid,GO,Evidence,Source)]
  
  check_dir(file.uniprot.go.selected)
  write.table(uniprot.filtered, file=file.uniprot.go.selected,quote=F,sep = "\t",row.names = F)

  myMessage("Uniprot database filtered.")
  
  setwd(currDir)
}


download_uniprot = function (id, trSize = 1000,keepnames=FALSE) {
  # id<-new_IDs[1:1000]
  # trSize<-1000
  if (!require(httr,quietly = T)) install.packages("httr")
  library(httr)
  currDir<-getwd()
  setwd(funDir)
  
  q<-0
  Nq<-length(id)
  if (file.exists(file.uniprot.download)) file.remove(file.uniprot.download)
  file.temp<-"temp.txt"
  while(q<Nq){
    items<-(q+1):min(q+trSize,Nq)
    write.table(id[items],file.temp,col.names = F,row.names = F,quote = F)
    
    r <- POST("https://www.uniprot.org/uploadlists/",
              body=
                list(file=upload_file(file.temp),
                     format= "txt",
                     from = "ACC+ID",
                     to = "ACC"),encode="multipart")
    text <- content(r, "text")
    check_dir(file.uniprot.download)
    write.table(text, file.uniprot.download,col.names = F,row.names = F,quote = F,append = T)
    q<-min(q+trSize,Nq)
    if (q %% 100==0 | q==Nq) message(paste0(Sys.time(),"  Genbank downloaded (",q,"/",Nq,")."))
    
    }
  
  fasta<-read_genbank(file.uniprot.download,output = "fasta",keepnames=keepnames)
  if (file.exists(file.temp)) file.remove(file.temp)
  
  setwd(currDir)
  return(fasta)
}

update_blastdb<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  if (!file.exists(file.uniprot.go.selected)) update_plantID()
  
  myMessage("Loading uniprot annotations.")
  uniprot.go<-fread(file.uniprot.go.selected)
  myMessage("Uniprot annotations loaded.")
  
  if (!file.exists(file.uniprot.fasta)){
    if (!file.exists(file.sprot.fasta)) download_sprot()
    myMessage("Loading Swissprot fasta.")
    if (!exists("read_fasta")) source("CommonUse.R")
    fasta.current<-read_fasta(file.sprot.fasta)
    myMessage("Swissprot fasta loaded.")
  }else{
    fasta.current<-read_fasta(file.uniprot.fasta)
    myMessage("Current fasta loaded.")
  }
  
  fasta.current[,ID:=sub(".*\\|(.*)\\|.*","\\1",fasta.current$ID)]
  fasta.current[,Text:=paste0(">",ID,"\n",Sequence)]
 
  new_IDs<-setdiff(uniprot.go$ID,fasta.current$ID)
  
  # fasta.noiea<-dowload_fasta(noiea_IDs[1:100])
  if (length(new_IDs)>0) {
    myMessage("Downloading neccessary fasta.")
    fasta.new<-download_uniprot(new_IDs,1000)
    myMessage("TrEMBL fasta downloaded.")
    fasta.merged<-rbind(fasta.current[ID %in% uniprot.go$ID],fasta.new)
    
    }
  else fasta.merged<-unique(fasta.current[ID %in% uniprot.go$ID])
  
  fasta.merged<-fasta.merged[!duplicated(ID)]
  
  check_dir(file.uniprot.fasta)
  write(fasta.merged$Text,file.uniprot.fasta)
  
  uniprot.nofasta<-unique(uniprot.go[!ID %in% fasta.merged$ID,.(ID,Taxid)])
  check_dir(file.uniprot.nofasta)
  write.table(uniprot.nofasta,file.uniprot.nofasta,sep = "\t",row.names = F,col.names = T,quote = F)
  
  myMessage("Blast database updated.")
  
  setwd(currDir)
}

prepare_blast1<-function(ver="5.1_with_Female"){
  if (!ver %in% c("3.1","Tak1 ver5.1","5.1_with_Female")) stop('Version should be "3.1","Tak1 ver5.1" or "5.1_with_Female".')
  if (!require(zip,quietly = T)) install.packages("zip")
  library(zip)
  currDir<-getwd()
  setwd(funDir)
  ver<-"5.1_with_Female"
  if (!file.exists(file.uniprot.fasta)) update_blastdb()
  if (ver=="3.1") {
    file.fasta.query<-file.fasta.ver3
  } else if (ver=="Tak1 ver5.1") {
      file.fasta.query<-file.fasta.ver5
  }else if (ver=="5.1_with_Female") {
        file.fasta.query<-file.fasta.ver5f}
  
  if (!file.exists(file.fasta.query)) {
    if (ver=="5.1_with_Female") merge_genome_ver5_with_female()
    else download_genome()
    }

  if (dir.exists(path.Blast1)) unlink(path.Blast1,recursive = T)
  check_dir(file.Blast1.query)
  file.copy(file.fasta.query,file.Blast1.query)
  
  file.copy(file.uniprot.fasta,path.Blast1)
  file.copy("Blast1.sh",path.Blast1)
  
  setwd(paste0(path.Blast1,".."))
  ziplist<-c(file.Blast1.query,file.uniprot.fasta,"Blast1.sh")
  # ziplist<-dir(path.Blast1, full.names = T)
  ziplist<-paste0(path.Blast1,sub(".*/(.*)$","\\1",ziplist))
  check_dir(file.Blast1.zip)
  zipr(file.Blast1.zip,ziplist,include_directories = T)
  message(paste0("Please use ",file.Blast1.zip," to do first blast search."))
  
  setwd(currDir)
}

prepare_blast2<-function(){
  if (!require(zip,quietly = T)) install.packages("zip")
  library(zip)
  currDir<-getwd()
  setwd(funDir)

  {
    # source("CommonUse.R")
  Hit_lim<-20000
  if (!file.exists(file.Blast1.result.self)) stop("Blast results (query-self) not found.")
  if (!file.exists(file.Blast1.result.uniprot)) stop("Blast results (query-uniprot) not found.")
  
  myMessage("Loading blast results")
  result.blast0<-read_blast(file.Blast1.result.self)
  result.blast1<-read_blast(file.Blast1.result.uniprot,cutoff=50)
  result.blast1[,qgene:=sub("\\..*","",query)]
  
  # file.uniprot<-"../../Blast2GO_Data/01-Parsed/uniprot_go_filtered.tsv"
  if (!file.exists(file.uniprot.go.selected)) update_plantID()
  myMessage("Loading uniprot database info.")
  info.uniprot<-fread(file.uniprot.go.selected)
  IDs.uniprot<-unique(info.uniprot[,.(ID,Taxid)])
  # result.blast1[,Taxid:=IDs.uniprot$Taxid[match(ID,IDs.uniprot$ID)]]
  if (sum(duplicated(IDs.uniprot$ID))>0) warning(paste("Multiple taxids for the same uniprot ID. Please check",file.uniprot))
  result.blast1<-merge(result.blast1,IDs.uniprot,all.x = T)
  
  if (!file.exists(file.fasta.query)) download_genome()
  myMessage("Loading query genome.")
  fasta.query<-read_fasta(file.fasta.query)
  IDs.query<-fasta.query[,.(Num=Num,gene=sub("\\..*","",ID),transcript=ID)]
  
  
  #Check clusters in query self-blast
  result.blast.clusters<-blast_cluster(result.blast0)
  blast.clusters<-data.table(separate_rows(result.blast.clusters,Genes,sep = ";"))
  blast.clusters<-merge(blast.clusters,IDs.query[,.(Genes=gene,query=transcript)],allow.cartesian = T)
  blast.clusters<-merge(blast.clusters,unique(result.blast1[,.(query,Hit=ID,Taxid)]),by="query",all.x=T,allow.cartesian = T)
  
  blast.clusters.sum<-blast.clusters[,.(Num_gene=length(unique(Genes)),
                                          Genes=paste(unique(Genes),collapse = ";"),
                                          Num_hit=length(unique(Hit[!is.na(Hit)])),
                                          Num_tax=length(unique(Taxid[!is.na(Taxid)]))
                                          ),by=Group]
  
  setorder(blast.clusters.sum,-Num_hit,-Num_tax,-Num_gene)
  
  # blast.clusters.renamed<-blast.clusters[Group %in% blast.clusters.sum[Num_hit>0]$Group]
  if (sum(blast.clusters.sum$Num_hit>Hit_lim)>0) warning("Size of single group exceeds limitation.")
  blast.clusters.sum[,NewGroup:=ifelse(Num_hit>Hit_lim,-1,0)]
  blast.clusters.sum[Num_hit==0,NewGroup:=-2]
  blast.clusters.sum[Num_hit>Hit_lim,NewGroup:=as.integer(factor(Group,levels = rev(Group)))]
  
    # blast.clusters.renamed[,Group:=as.integer(factor(Group,levels=unique(Group)))]
  myMessage("Merging small clusters.")
  grList<-unique(blast.clusters.sum[NewGroup==0,Group])
  cluster<-grList[1]
  gct<-max(blast.clusters.sum$NewGroup)+1
  for (g in 1:length(grList)){
    cluster<-union(cluster,grList[g])
    Nhit<-length(unique(blast.clusters[!is.na(Hit) & Group %in% cluster]$Hit))
    if (g==length(grList)) 
      blast.clusters.sum[Group %in% cluster, NewGroup:=gct] 
    else if (Nhit>Hit_lim){
      blast.clusters.sum[Group %in% setdiff(cluster,g), NewGroup:=gct] 
      gct<-gct+1
      cluster<-grList[g]
    }
  }
  
  
  blast.clusters.new<-merge(blast.clusters,blast.clusters.sum[NewGroup>0,.(Group,NewGroup)],by="Group")
  blast.clusters.new[,`:=`(Group=NewGroup,
                           NewGroup=NULL
                           )]
  }
  blast.clusters.ids<-unique(rbind(blast.clusters.new[!is.na(Hit),.(Group,ID=Hit,Taxid)],blast.clusters.new[,.(Group,ID=query,Taxid=0)]))
  setorder(blast.clusters.ids,Group,Taxid,ID)
  blast.clusters.ids[,Tax_idx:=as.integer(factor(Taxid,levels=unique(Taxid)))-1,by=Group]
  blast.clusters.ids[,Seq_idx:=1:.N,by=c("Group","Tax_idx")]

  #Split uniprot fasta by group
  if (!file.exists(file.uniprot.fasta)) update_blastdb()
  myMessage("Loading uniprot fasta.")
  fasta.uniprot<-read_fasta(file.uniprot.fasta)
  
  check_dir(file.Blast2.hits)
  fasta.hits<-fasta.uniprot[ID %in% blast.clusters.ids$ID]
  write(unique(fasta.hits$Text),file.Blast2.hits)
  
  check_dir(file.Blast2.groups)
  write.table(blast.clusters.ids,file.Blast2.groups,quote = F,row.names = F)
  
  for (gr in unique(blast.clusters.ids$Group)){
    # gr<-1
    group.ids<-blast.clusters.ids[Group==gr]
    group.uniprot<-fasta.uniprot[ID %in% group.ids$ID]
    group.uniprot<-merge(group.uniprot,group.ids[,.(ID,Tax_idx,Seq_idx)],by="ID",all.x = T)
    setorder(group.uniprot,Tax_idx,Seq_idx)
    group.uniprot[,Idx:=paste0(Tax_idx,"_",Seq_idx)]
    
    group.uniprot[,NewText:=paste0(">",Idx,"\n",Sequence)]
    group.uniprot.fasta<-paste0(path.Blast2,"Group_",gr,".fasta")
    # print(setdiff(group.uniprot$ID,group.ids$ID))
    write(unique(group.uniprot$NewText),group.uniprot.fasta)
  }
  
  file.copy("Blast2.sh",path.Blast2)
  file.copy(file.fasta.query,file.Blast2.query)
  # file.copy(file.uniprot.fasta,Blast2)
  
  setwd(paste0(path.Blast2,".."))
  ziplist <- dir(path.Blast2, full.names = TRUE)
  ziplist<-setdiff(gsub("/+","/",ziplist),file.Blast2.zip)
  check_dir(file.Blast2.zip)
  zipr(file.Blast2.zip,ziplist,include_directories = T)
  
  #Splitting blast1
  # result.blast1.gr<-merge(result.blast1,blast.clusters.new[,.(query,ID=Hit,Group)],by=c("query","ID"),all.x=T,allow.cartesian = T)
  

  setwd(currDir)
}

prepare_orthofinder<-function(){
  {
  currDir<-getwd()
  setwd(funDir)

  check_dir(path.ortho)
  if (!file.exists(file.Blast2.groups)) prepare_blast2()
  
  IDs<-fread(file.Blast2.groups)
  IDs[,Idx:=paste0(Tax_idx,"_",Seq_idx)]
  Groups<-unique(IDs$Group)
  blast.q2q<-read_blast(file.Blast1.result.self)
  blast.q2u<-read_blast(file.Blast1.result.uniprot)
  blast.u2q<-read_blast(file.Blast2.result.query)
  blast.oriID<-rbind(rbind(blast.q2q,blast.q2u),blast.u2q)
  header.blast<-c("query", "ID", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  fasta.uniprot<-read_fasta(file.uniprot.fasta)
  fasta.query<-read_fasta(file.fasta.query)
  fasta.all<-rbind(fasta.query,fasta.uniprot)
  
  
  for (gr in Groups){
    path.gr<-paste0(path.ortho,"Group_",gr,"/")
    check_dir(path.gr)
    IDs.gr<-IDs[Group==gr]
    file.seqids<-paste0(path.gr,"SequenceIDs.txt")
    write(paste0(IDs.gr$Idx,": ",IDs.gr$ID),file.seqids)
    file.speids<-paste0(path.gr,"SpeciesID.txt")
    write(unique(paste0(IDs.gr$Tax_idx,": ",IDs.gr$Taxid,".fasta")),file.speids)
    blast.oriID.gr<-blast.oriID[query %in% IDs.gr$ID & ID %in% IDs.gr$ID]
    blast.oriID.gr[,`:=`(q.idx=IDs.gr$Idx[match(query,IDs.gr$ID)],
                         h.idx=IDs.gr$Idx[match(ID,IDs.gr$ID)]
                         )]
    blast.qu.gr<-blast.oriID.gr[,.(query=q.idx,
                                   ID=h.idx,
                                   pident, length, mismatch, gapopen, qstart, 
                                   qend, sstart,send, evalue,bitscore
                                   )]
    file.blast.gr<-paste0(path.Blast2.results,"Blast-Group_",gr,".fasta.tsv")
    blast.self.gr<-read_blast(file.blast.gr)
    blast.gr<-rbind(blast.qu.gr,blast.self.gr[,..header.blast])
    blast.gr[,`:=`(qtax=sub("_.*","",query),
                   htax=sub("_.*","",ID)
                   )]
    fasta.gr<-fasta.all[ID %in% IDs.gr$ID]
    fasta.gr<-merge(fasta.gr,IDs.gr,by="ID",all.x=T,allow.cartesian = T)
    fasta.gr[,NewText:=paste0(">",Idx,"\n",Sequence)]
    for (qt in unique(IDs.gr$Tax_idx)){
      for (ht in unique(IDs.gr$Tax_idx)){
        file.out<-paste0(path.gr,"Blast",qt,"_",ht,".txt")
        
        blast.sub<-blast.gr[qtax==qt & htax==ht]
        write.table(blast.sub[,..header.blast],file.out,quote = F,col.names = F,row.names = F,sep = "\t")
      }
      fasta.qt<-unique(fasta.gr[Tax_idx==qt])
      file.out.fasta<-paste0(path.gr,"Species",qt,".fa")
      write(fasta.qt$NewText,file.out.fasta)
    }
    myMessage(paste0("File splitting finished for Group ",gr))
  }
  }
  setwd(currDir)
}
parse_blast_by_species<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  file.fasta<-"../../Blast2GO_Data/01-Parsed/Blast/uniprot_selected.fasta"
  if (!file.exists(file.fasta)) update_blast()
  
  myMessage("Loading uniprot fasta file.")
  fasta<-read_fasta(file.fasta)
  fasta[,`:=`(ID=sub(".*\\|(.*)\\|.*","\\1",ID),
              Taxid=sub(" .*","",sub(".*OX=","",Text)),
              Taxon=sub(" OX.*","",sub(".*OS=","",Text))
              )]
  
  file.blast<-"../../Blast2GO_Data/01-Parsed/Blast_Mpa_v3.1-uniprot.tsv"
  if (!file.exists(file.blast)) stop("Blast file not found.")
  
  myMessage("Loading blast results.")
  blast<-fread(file.blast)
  names(blast)<-c("query", "ID", "pident", "ppos", "evalue", "bitscore")
  setorder(blast,query,evalue,-bitscore,-ppos)
  blast[,hit_num:=1:.N,by=query]
  blast<-merge(blast,unique(fasta[,.(ID,Taxid)]),all.x = T)
  blast.tax<-blast[,.(Num=.N,IDs=paste0(ID,collapse = ";")),by=.(query,Taxid)]
  blast.tax.1<-blast[,.(taxid_num=length(unique(Taxid))),by=query]
  
  fasta.sub<-fasta[ID %in% blast.tax[Num>1,ID]]
  fasta.stats<-fasta.sub[,.(Num=.N),by=.(Taxid,Taxon)]
  setorder(fasta.stats,-Num,Taxon)
  fasta.stats[,AccPer:=cumsum(Num)/sum(Num)]
  file.stats<-"../../Blast2GO_Data/01-Parsed/Blast/uniprot_selected_species_dist.tsv"
  write.table(fasta.stats,file.stats,row.names = F,col.names = T,sep="\t")
  
  test<-blast[ID %in% fasta$ID[fasta$Taxid %in% fasta.stats$Taxid[fasta.stats$Num==2]]]
  
  fastaPATH<-"../../Blast2GO_Data/01-Parsed/Blast/BySpecies"
  if (dir.exists(fastaPATH)) unlink(fastaPATH,recursive = T)
  dir.create(fastaPATH)
  for (s in fasta.stats[Num>1,Taxid]){
    fasta.trunk<-fasta.sub[Taxid==s,Text]
    file.trunk<-paste0(fastaPATH,"/",s,".fasta")
    write.table(fasta.trunk,file.trunk,row.names = F,col.names = F)
  }
  myMessage("Blast database splitted by species.")
  file.Mp<-"../../Blast2GO_Data/00-Downloaded/Mpolymorphav3.1.allTrs.pep_annot.fa"
  if (!file.exists(file.Mp)) download_genome()
  
  file.copy(file.Mp,paste0(fastaPATH,"/Mp.fasta"))
  setwd(currDir)
}

my_blast<-function(){
  currDir<-getwd()
  setwd(funDir)
  sysType<-Sys.info()[['sysname']]
  if (sysType=="Windows") {
    nCPU<-as.numeric(system("WMIC CPU Get NumberOfCores",intern = T)[2])
    system("setx BLASTDB_LMDB_MAP_SIZE 1000000")
  }
  if (sysType=="Darwin") nCPU<-as.numeric(system("getconf _NPROCESSORS_ONLN",intern = T))
  
  
  dbPATH<-"../../Blast2GO_Data/01-Parsed/Blast/"
  setwd(dbPATH)
  file.db<-"uniprot_selected.fasta"
  system(paste0("makeblastdb -in ",file.db," -parse_seqids -dbtype prot -out ", file.db))
  setwd("../")
  
  myMessage("Blast database built.")
  
  file.query<-"../00-Downloaded/Mpolymorphav3.1.allTrs.pep_annot.fa"
  file.test<-"GA+SABBATH.fasta"
  blastcmd<-paste0( 'blastp -db ./Blast/',file.db,
                   ' -query ',file.query,
                   # ' -query ',file.test,
                   ' -out Blast_Mpa_v3.1-uniprot.tsv -evalue 1e-5',
                   ' -outfmt "6 qaccver saccver pident ppos evalue bitscore"',
                   ' -num_threads ',nCPU)
  system(blastcmd)
  
  setwd(currDir)
}


