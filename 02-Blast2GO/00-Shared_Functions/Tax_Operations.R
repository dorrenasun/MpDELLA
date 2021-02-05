if (!require(data.table,quietly =T)) install.packages("data.table")
library(data.table)

NothingButToDetectDirectory<-function(){}
funDir<-paste0(getwd(),"/",getSrcDirectory(NothingButToDetectDirectory))
# funDir<-getSrcDirectory(NothingButToDetectDirectory)

if (!exists("myMessage")) source(paste0(funDir,"/CommonUse.R"))
source(paste0(funDir,"/File_Locations.R"))

download_taxonomy<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading NCBI taxonomy.")
  check_dir(file.taxdmp)
  download.file(url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip",
                # method = "wget",
                destfile = file.taxdmp)
  
  message(paste0("Taxonomy downloaded at: ",file.info(file.taxdmp)$ctime))
  myMessage("NCBI taxonomy downloaded.")
  setwd(currDir)
}


update_taxonomy<-function(){
  currDir<-getwd()
  setwd(funDir)
 { 
  myMessage("Updating taxonomy.")
  if (!file.exists(file.taxdmp)) download_taxonomy()
  taxPATH<-"../../Blast2GO_Data/01-Parsed/Taxonmy/"
  if (!dir.exists(taxPATH)) dir.create(taxPATH,recursive = T)
  unzip(file.taxdmp,
        files = c("nodes.dmp","names.dmp","merged.dmp"),
        exdir = taxPATH
  )
  
  file.nodes<-paste0(taxPATH,"nodes.dmp")
  nodes<-fread(file.nodes,header = F,sep = "|")
  nodes[,`:=`(V14=NULL,
              V3=gsub("\t","",V3,fixed = T),
              V4=gsub("\t","",V4,fixed = T)
  )]
  
  names(nodes)<-c("taxid", "parent", "rank", "embl_code", "division_id", 
                  "div_flag", "genetic_code", "GC_flag", "mt_genetic_code",
                  "MGC_flag", "GenBank_flag", "subtree_flag", "comments")
  
  file.names<-paste0(taxPATH,"names.dmp")
  names<-fread(file.names,header = F,sep = "|")
  names[,`:=`(V5=NULL,
              V2=gsub("\t","",V2,fixed = T),
              V3=gsub("\t","",V3,fixed = T),
              V4=gsub("\t","",V4,fixed = T))]
  names(names)<-c("taxid", "name", "unique_name", "class")
  
  file.merged<-paste0(taxPATH,"merged.dmp")
  merged<-fread(file.merged,header = F,sep = "|")
  merged[,V3:=NULL]
  names(merged)<-c("oldid","taxid")
  
  taxon<-merge(nodes[,.(taxid,parent,rank)],names[class=="scientific name",.(taxid,name)],all.x = T)
  taxon<-merge(merged,taxon,by="taxid",all = T)
}
  taxon[taxid==1,rank:="root"]

  check_dir(file.taxon)
  write.table(taxon, file=file.taxon,quote=F,sep = "\t",row.names = F,na = "")
  
  file.remove(c(file.nodes,file.names,file.merged))
  setwd(currDir)
  myMessage("Taxonomy updated.")
}

load_taxonomy<-function(){
  currDir<-getwd()
  setwd(funDir)
  #Load taxon from file
  if (!file.exists(file.taxon)) update_taxonomy()
  myTaxon<-fread(file.taxon)
  myTaxon[taxid==1, rank:="root"]
  assign("myTaxon",myTaxon,envir = .GlobalEnv)
  myMessage("NCBI taxonomy loaded.")
  setwd(currDir)
}

get_tax_descendants<-function(query){
  currDir<-getwd()
  setwd(funDir)
  if (!file.exists(file.taxon)) update_taxonomy()
  
  myTaxon<-fread(file.taxon)
  query<-as.numeric(query)
  myMessage(paste0("Getting taxonomy descendants for: ",myTaxon[taxid==query,name]))
  
  myParent<-query
  n<-0
  while (n!=length(myParent)){
    n<-length(myParent)
    child<-myTaxon[parent %in% myParent,taxid]
    myParent<-union(myParent,child)
  }
  myTaxon.sub<-myTaxon[taxid %in% myParent]
  myMessage("Taxonomy descendants get.")
  
  setwd(currDir)
  return(myTaxon.sub)
}

get_suptax<-function(rank){
  if (!exists("myTaxon")) load_taxonomy()
  if (!exists("TaxRanks")){
    taxon<-myTaxon[,.(taxid,parent,rank)]
    taxon[,prank:=rank[match(parent,taxid)]]
    TaxRanks<-unique(taxon[!(rank=="no rank" | prank=="no rank" | rank==prank),.(rank,prank)])
    assign("TaxRanks",TaxRanks,envir = .GlobalEnv)
  }
  des<-rank
  n<-0
  while (n!=length(des)){
    n<-length(des)
    parent<-TaxRanks[rank %in% des,prank]
    des<-union(des,parent)
  }
  return(unique(des))
}

fix_tax_IDs<-function(taxids,sum_at=NULL){
  currDir<-getwd()
  setwd(funDir)
  if (!exists("myTaxon")) load_taxonomy()
  # sum_at<-"species"
  if (!is.null(sum_at)){
    if (!sum_at %in% myTaxon$rank){
      message("'sum_at' should be a legal taxonomy rank. Please choose one of the following:")
      message(paste(unique(myTaxon$rank),collapse = ", "))
      stop("Illegal 'sum_at' argument.")
    }
  }
  
  # taxids<-uniprot.all$Taxid
  dt<-data.table(Taxid=taxids)
  dt[,updated:=ifelse(Taxid %in% myTaxon$oldid,myTaxon$taxid[match(Taxid,myTaxon$oldid)],Taxid)]
  
  dt[,rank:=myTaxon[match(updated,taxid),rank]]
  
  if (is.null(sum_at)) {
    setwd(currDir)
    return(dt$updated)
  }else {
    dt[,`:=`(sum=updated,srank=rank)]
    
    suptax<-get_suptax(sum_at)
    while (nrow(dt[(!srank %in% suptax) & updated %in% myTaxon$taxid])>0){
      dt[(!srank %in% suptax) & updated %in% myTaxon$taxid,sum:=myTaxon[match(sum,taxid),parent]]
      dt[(!srank %in% suptax) & updated %in% myTaxon$taxid,srank:=myTaxon[match(sum,taxid),rank]]
    }
    setwd(currDir)
    return(dt$sum)
  }
  
}
