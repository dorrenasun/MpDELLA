if (!require(data.table,quietly =T)) install.packages("data.table")
library(data.table)
if (!require(tidyr,quietly = T)) install.packages("tidyr")
library(tidyr)
if (!require(R.utils,quietly = T)) install.packages('R.utils')

myMessage<-function(str,quiet=FALSE){
  # quiet<-TRUE
  quiet<-ifelse(exists("quiet_default"),quiet_default,quiet)
  if (!require(beepr)) install.packages("beepr")
  library(beepr)
  message(paste0(Sys.time(),"  ",str))
  if (!quiet==TRUE) {
    if (Sys.info()[['sysname']]=="Darwin") {
      system(paste0("say ",str))
    }
    else beep(1)
  }
}

check_dir<-function(filename){
  path<-sub("(.*)\\/.*","\\1",filename)
  if (!dir.exists(path)) dir.create(path,recursive = T)
}

read_fasta<-function(file){
  library(R.utils)
  
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

read_genbank<-function(file,output="fasta",keepnames=FALSE){
  if (!file.exists(file)) stop("Genbank file not found.")
  if (!output %in% c("fasta","GO")) stop('Output should be "fasta" or "GO".')
  # file<-"downloaded.db"
  myMessage("Loading genbank file...")
  # file<-file.uniprot.download
  gb.ori<-fread(file,fill=T,sep=" ",header = F,strip.white = F)
  
  
  myMessage("Selecting tags...")
  gb.ori[,Tag:=sub(" +.*","",V1)]
  gb.ori[Tag=="",Tag:="SEQ"]
  Tag.select<-c("ID","AC","DE","OX","DR","SQ","SEQ")
  gb.trans<-gb.ori[Tag %in% Tag.select]
  gb.trans[,Content:=substring(V1,6,length(V1))]
  
  myMessage("Fetching IDs...")
  ID.positions<-which(gb.trans$Tag=="ID")
  ID.lengths<-c(ID.positions[-1],nrow(gb.trans)+1)-ID.positions
  gb.trans[,Num:=rep(1:length(ID.positions),ID.lengths)]
  gb.trans<-dcast(gb.trans, Num ~ Tag, fun.aggregate = function(x){paste0(x,collapse = "!")},value.var = "Content")
  
  myMessage("Formatting...")
  gb.trans<-gb.trans[,.(Num=Num,
                         ID=sub(";.*","",AC),
                         Taxid=sub(".*=([0-9]+)[; ].*","\\1",OX),
                        UniprotID=sub(" .*","",ID),
                         DE=sub("RecName: Full=","",sub("[;!{}].*","",DE)),
                         DR=DR,
                         SEQ=gsub("!","\n",SEQ,fixed = T),
                         Length=as.integer(sub(".* ([0-9]+) .*","\\1",sub(";.*","",SQ)))
                                   
  )]
  if (sum(!nchar(gsub("[ |\\\n]","",gb.trans$SEQ))==gb.trans$Length)>0) warning("Something went wrong in parsing sequences.")
  if (output=="fasta"){
    myMessage("Extracting fasta.")
    if (keepnames==T){
      gb.fasta<-gb.trans[,.(Num=Num,
                            ID=ID,
                            UniprotID=UniprotID,
                            DE=DE,
                            Sequence=gsub("[ |\\\n]","",SEQ),
                            Text=paste0(">",ID,"\n",SEQ)
                            )]

    }else{
      gb.fasta<-gb.trans[,.(Num=Num,
                            ID=ID,
                            Sequence=gsub("[ |\\\n]","",SEQ),
                            Text=paste0(">",ID,"\n",SEQ)
                            )]
      
    }
    return(gb.fasta)
  }
  
  if (output=="GO"){
    myMessage("Extracting GOs.")
    gb.trans<-data.table(separate_rows(gb.trans,DR,sep = "!"))
    gb.trans[,DR_Tag:=sub(";.*","",DR)]
    
    gb.go<-gb.trans[DR_Tag=="GO",.(Num,ID,Taxid,DR,DR_Tag)]
  
    gb.go<-gb.go[,c("Tag","GO","Name","Evidence"):=tstrsplit(DR,";")]
    gb.go[,c("Evidence","Source"):=tstrsplit(Evidence,"[:\\.]")]
    gb.go[,GO:=sub(" ","",GO)]
    return(gb.go[,.(ID,Taxid,GO,Evidence,Source)])
  }
}

read_gaf<-function(file){
  if (!file.exists(file)) stop("Gaf file not found.")
  # file<-file.noiea.gaf
  gaf.ori<-fread(file,header=F,strip.white = F)
  header.gaf<-c("DB", "DB_ID", "DB_Symbol", 
                "Qualifier", "GO", "DB_Ref", "Evidence", 
                "With_From", "Aspect", "DB_Name", 
                "Synonym", "DB_Type", "Taxon", "Date", 
                "Source", "Extension", "Gene_Product_ID")
  names(gaf.ori)<-header.gaf
  # gaf.ori[,Accession:=ifelse(DB_ID %in% sprot.go$ID,sprot.go$Accession[match(DB_ID,sprot.go$ID)],DB_Symbol)]
  gaf.go<-gaf.ori[Qualifier=="",.(ID=DB_ID,
                                      Taxid=sub("taxon:([0-9]+).*","\\1",Taxon),
                                      GO=GO,
                                      Evidence=Evidence,
                                      Source=Source
  )]
  return(gaf.go)
}

read_blast<-function(file, cutoff=Inf,by="query"){
  if (!file.exists(file)) stop("Blast result file not found.")
  if (!by %in% c("query","gene")) stop("'by' should be either 'query' or 'gene'.")
  # file<-"../../Blast2GO_Data/01-Parsed/Blast_results/Blast_query-self.tsv"
  # cutoff<-30
  header.blast<-c("query", "ID", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "ppos")
  
  blast<-fread(file,header = F)
  names(blast)<-header.blast
  setorder(blast,query,-bitscore,evalue,-ppos,-pident)
  if (by=="gene") {
    blast[,hgene:=sub("\\..*","",ID)]
    blast[,hit_num:=as.numeric(factor(hgene,levels = unique(hgene))),by=query]
  }
  if (by=="query") blast[,hit_num:=1:.N,by=query]
  return(blast[hit_num<=cutoff])
}

blast_cluster<-function(blast_result){
  # blast_result<-result.blast0
  blast_result[,`:=`(qgene=sub("\\..*","",query),
                hgene=sub("\\..*","",ID))]
  blast_result.clean<-unique(blast_result[,.(qgene,hgene)])
  blast_result.sym<-unique(rbind(blast_result.clean,blast_result.clean[,.(qgene=hgene,hgene=qgene)]))
  # blast_result.clean<-data.table(qgene=c("B","C","C","C","D","D","F"),hgene=c("C","B","D","E","C","E","A"))
  
  left<-unique(blast_result.sym$qgene)
  total<-length(left)
  message(paste0("Total number of genes:", total))
  count<-1
  prog_ct<-0
  myMessage("Start finding clusters.")
  message(paste0("|",paste0(rep("-",50),collapse = ""),"|"))
  message("|",appendLF = F)
  while(length(left)>0){
    node<-left[1]
    n<-0
    while (n!=length(node)){
      n<-length(node)
      hit<-blast_result.sym[qgene %in% node,hgene]
      node<-union(node,hit)
    }
    temp<-data.table(Group=count,
                     Genes=paste0(node,collapse = ";"),
                     Num=length(node))
    if (count==1) {
      group<-temp
      } else group<-rbind(group,temp)
    left<-setdiff(left,node)
    count<-count+1
    prog<-(total - length(left)) %/% (total/50)
    if ((prog - prog_ct)>0) message(paste0(rep("=",prog - prog_ct),collapse = ""),appendLF = F)
    prog_ct<-prog
    # message(paste0("There are ",length(left)," genes left."))
  }
  message("|")
  setorder(group,-Num)
  group[,Group:=1:nrow(group)]
  # group<-data.table(separate_rows(group,Genes,sep = ";"))
  myMessage("Clusters identified from blast.")
  return(group)
}

seg_fasta<-function(file,trSize=5000,outpath=NULL){
  # file<-"../../Blast2GO_Data/00-Downloaded/Mpolymorphav3.1.allTrs.pep_annot.fa"
  # outpath<-NULL
  # trSize<-5000

  inFile.name.split<-unlist(strsplit(sub(".*\\/","",file),".",fixed = T))
  for (n in length(inFile.name.split):1){
    inFile.name.ext<-inFile.name.split[n]
    if (! inFile.name.ext %in% c("gz","bz2")) {
      inFile.name.base<-paste0(inFile.name.split[1:n-1],collapse = ".")
      inFile.name.ext<-paste0(".",inFile.name.ext)
      break()
    }
  }

  inFile.path<-sub("(.*)\\/.*","\\1",file)
  if (!file.exists(file)) stop("Fasta file not found!")
  if (is.null(outpath)) outpath<-paste0(inFile.path,"/Fasta_segs/")
  if (! substr(outpath,nchar(outpath),nchar(outpath)) %in% c("/","\\")) outpath<-paste0(outpath,"/")
  if (dir.exists(outpath)) {
    unlink(substr(outpath,1,nchar(outpath)-1),recursive = T)
  }
  dir.create(outpath)
  
  data<-read_fasta(file)
  data.unique<-data[!duplicated(Sequence)]
  data.duplicated<-merge(data[duplicated(Sequence),.(ID,Sequence)],data.unique[,.(ID,Sequence)],by="Sequence",all.x = T)
  names(data.duplicated)<-c("Sequence","ID","DuplicateOf")
  outFile.duplicate<-paste0(outpath,"Duplicate_Sequences.tsv")
  write.table(data.duplicated,outFile.duplicate,row.names = F,quote = F,sep = "\t")
  
  for (i in 1:ceiling(nrow(data.unique)/trSize)){
    outFile<-paste0(outpath,inFile.name.base,"_",i,inFile.name.ext)
    items<-(trSize*(i-1)+1):min(trSize*i,nrow(data.unique))
    # print(paste(outFile,"\n",min(items),max(items),length(items),sep = ", "))
    write.table(data.unique[items,Text],outFile,row.names = F,col.names = F,quote = F)
  }
  message(paste0("Trunctated fasta files saved to: ",outpath))
}

