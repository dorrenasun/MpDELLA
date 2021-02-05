if (!require(ontologyIndex,quietly =T)) install.packages("ontologyIndex")
library(ontologyIndex)
if (!require(data.table,quietly =T)) install.packages("data.table")
library(data.table)
if (!require(tidyr,quietly =T)) install.packages("tidyr")
library(tidyr)
if (!require(beepr,quietly =T)) install.packages("beepr")

file.GO<-"./go.obo"

# currDir<-getwd()
NothingButToDetectDirectory<-function(){}
funDir<-paste0(getwd(),"/",getSrcDirectory(NothingButToDetectDirectory))

check_dir<-function(filename){
  path<-sub("(.*)\\/.*","\\1",filename)
  if (!dir.exists(path)) dir.create(path,recursive = T)
}

download_GO<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  myMessage("Downloading gene ontologies.")
  check_dir(file.GO)
  download.file(url = "http://purl.obolibrary.org/obo/go.obo",destfile = file.GO)
  myMessage("Gene ontologies downloaded.")
  
  setwd(currDir)
}

load_go<-function(){
  currDir<-getwd()
  setwd(funDir)
  
  file.OBO<-file.GO
  if (!file.exists(file.OBO)) download_GO()
  
  #Load go.obo
  list.relations<-get_relation_names(file.OBO)
  list.relations<-c("is_a","part_of","regulates","negatively_regulates","positively_regulates","occurs_in")
  myOntology<-get_ontology(file = file.OBO,propagate_relationships = list.relations,extract_tags = "everything")
  assign("myOntology",myOntology,envir = .GlobalEnv)
  
  myMessage("Gene ontologies loaded.")
  setwd(currDir)
  
}


fix_alternative<-function(dt,fix="GO"){
  
  #Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  #Find all GOs with alternative IDs
  if (!exists("GO.list.alt")){
    GO.list.alt<-data.table(GO=myOntology$id,Alt=lapply(myOntology$alt_id, function(x)paste0(x,collapse = ",")))
    GO.list.alt<-data.table(separate_rows(GO.list.alt[!Alt==""],Alt,sep = ","))
    assign("GO.list.alt",GO.list.alt,envir = .GlobalEnv)
  }

  # Check for legal input
  if (sum(fix %in% names(dt))<1) {
    msg<-paste0("Change the 'fix' argument or check the input colomn names for '",fix,"'.",collapse = "")
    stop(msg)
  }
  
  # Fix alternative IDs
  if( sum(!unlist(dt[,..fix]) %in% myOntology$id)==0) {
    message("No GOs in alternative ID detected!")
  } else{
    dt$GO_blabla<-dt[,..fix]
    
    message(paste0("GO in alternative IDs: ",sum(!dt$GO_blabla %in% myOntology$id)))
    dt[GO_blabla %in% GO.list.alt$Alt,Fixed:=GO.list.alt$GO[match(GO_blabla,GO.list.alt$Alt)]]
    dt[is.na(Fixed),Fixed:=GO_blabla]
    message(paste0("After fixation: ",sum(!dt$Fixed %in% myOntology$id)))
    
    OriginalName<-"OriginalGO"
    ct<-1
    while (OriginalName %in% names(dt)){
      OriginalName<-paste0("OriginalGO","_",ct,collapse = "")
      ct<-ct+1
    }
    setnames(dt,fix,OriginalName)
    setnames(dt,"Fixed",fix)
    dt[,"GO_blabla":=NULL]
  }
  myMessage("Alternative IDs fixed!")
  return(dt)
}

  
fix_obsolete<-function(dt,fix="GO"){
  dt<-list.sub
  fix<-"GO"
  #Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  #Find all obsolete GOs
  if (!exists("GO.list.obs")){
    GO.list.obs<-data.table(GO=names(myOntology$obsolete)[myOntology$obsolete])
    GO.list.obs[,`:=`(Ont=get_prop(GO,"namespace"),
                      Name=get_prop(GO,"name"),
                      Replace=as.character(get_prop(GO,"replaced_by")))]
    GO.list.obs[Replace=="character(0)",Replace:=""]
    assign("GO.list.obs",GO.list.obs,envir = .GlobalEnv)
  }

  # Check for legal input
  if (sum(fix %in% names(dt))<1) {
    msg<-paste0("Change the 'fix' argument or check the input colomn names for '",fix,"'.",collapse = "")
    stop(msg)
  }
  
  # Remove obsolete terms
  dt$GO_blabla<-dt[,..fix]
  
  NumObs<-sum(unique(dt[!is.na(GO_blabla),GO_blabla]) %in% GO.list.obs$GO)
  if( NumObs==0) {
    dt[,GO_blabla:=NULL]
    message("No obsolete GO detected!")
  } else{
    message(paste0("Obsolete GOs: ",NumObs))
    
    dt[GO_blabla %in% GO.list.obs$GO,Replace:=GO.list.obs$Replace[match(GO_blabla,GO.list.obs$GO)]]
    
    list.replace<-unique(dt[!Replace=="",.(GO_blabla,Replace)])
    list.replace[,msg:=paste0(GO_blabla," to ",Replace)]
    message(paste0("Obsolete GOs replaced:\n  ",paste0(list.replace$msg,collapse = "\n  ")))
    dt[is.na(Replace),Replace:=GO_blabla]
    
    list.remove<-unique(dt[Replace=="",GO_blabla])
    
    dt<-dt[!Replace==""]
    message(paste0("Obsolete GOs removed:\n  ",paste0(list.remove,collapse = ",  ")))
    
    dt[,c(fix,"GO_blabla"):=NULL]
    setnames(dt,"Replace",fix)
  }
  return(dt)
}

remove_redundancy<-function(dt,trSize=100,select.by="query",go="GO"){
  # Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  # Check for legal input
  if (sum(c(select.by,go) %in% names(dt))<length(select.by)+1) {
    msg<-paste0("Check the input colomn names for ",paste0("'",select.by,"'",collapse = ", ")," and 'GO'.",collapse = "")
    stop(msg)
  }
  
  # Remove redundancy
  myMessage("Start marking redundancy...")

  setkeyv(dt,select.by)
  query.list<-unique(dt[,..select.by])
  query.list[,unique_id_blabla:=1:nrow(query.list)]
  dt<-merge(dt,query.list,by=select.by)
  dt$GO_blabla<-dt[,..go]

  q<-0
  dt$Flag_blabla<-NA
  while(q<nrow(query.list)){
    items<-which(dt$unique_id_blabla %in% query.list[(q+1):min(q+trSize,nrow(query.list)),unique_id_blabla])
    dt[items, Flag_blabla:=GO_blabla %in% minimal_set(myOntology,GO_blabla),by=select.by]
    q<-min(q+trSize,nrow(query.list))
    if (q %% 100==0 | q==nrow(query.list)) message(paste0(Sys.time(),"  Redundancy marked (",q,"/",nrow(query.list),")."))
  }
  dt.out<-dt[Flag_blabla==TRUE]
  dt.out[,unique_id_blabla:=NULL]
  dt.out[,Flag_blabla:=NULL]
  dt.out[,GO_blabla:=NULL]
  if (exists("quiet_default")) {
    if (!quiet_default) beep(1)
  }else beep(1)
  
  myMessage("Redundancy removed!")
  return(dt.out)
}

find_anc<-function(golist){
  # Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  GO.list.ori<-data.table(GO=unique(golist))
  GO.list.ori[,Anc:=sapply(GO,function(x){paste0(get_ancestors(myOntology,x),collapse = ",")})]
  GO.list.anc<-data.table(separate_rows(GO.list.ori,Anc,sep = ","))
  GO.list.anc<-GO.list.anc[!Anc==""]
  # GO.list.anc<-GO.list.anc[Anc %in% GO]
  GO.list.anc$Tag<-""
  GO.list.anc[Anc==GO, Tag:="Ori"]
  
  return(GO.list.anc)
}

find_dec<-function(golist){
  # Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  GO.list.ori<-data.table(GO=unique(golist))
  GO.list.ori[,Dec:=sapply(GO,function(x){paste0(get_descendants(myOntology,x),collapse = ",")})]
  GO.list.dec<-data.table(separate_rows(GO.list.ori,Dec,sep = ","))
  GO.list.dec<-GO.list.dec[!Dec==""]
  # GO.list.anc<-GO.list.anc[Anc %in% GO]
  GO.list.dec$Tag<-""
  GO.list.dec[Dec==GO, Tag:="Ori"]
  
  return(GO.list.dec)
}


num_anc<-function(golist){
  # Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  GO.list.ori<-data.table(GO=unique(golist))
  GO.list.ori[,NumAnc:=sapply(GO,function(x)length(get_ancestors(myOntology,x)))]
  return(GO.list.ori[match(golist,GO),NumAnc])
}

get_prop<-function(golist,property){
  # Check if GOs were loaded
  if (!exists("myOntology")) load_go()
  
  GO.list<-data.table(GO=unique(golist))
  GO.list[,X:=sapply(GO,function(x) get_term_property(myOntology,property,x))]
  return(GO.list[match(golist,GO),X])
}



# setwd(currDir)

