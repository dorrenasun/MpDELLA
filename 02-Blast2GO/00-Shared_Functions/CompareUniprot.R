source("./UpdateDatabase.R")
source("./GO_Operations.R")


file.gaf<-"../../Blast2GO_Data/00-Downloaded/tair.gaf.gz"
gaf.At<-fread(file.gaf,header=F,strip.white = F)
gaf.At[,ID:=sub("\\|.*","",Synonym)]
names(gaf.At)<-c("DB",
                 "DB Object ID",
                 "DB Object Symbol",
                 "Qualifier",
                 "GO",
                 "DB:Reference",
                 "Evidence Code",
                 "With (or) From",
                 "Aspect",
                 "DB Object Name",
                 "Synonym",
                 "DB Object Type",
                 "Taxon",
                 "Date",
                 "Assigned By",
                 "Annotation Extension",
                 "Gene Product Form ID")
gaf.At.not<-gaf.At[Qualifier=="NOT"]

At_ids<-get_tax_descendants(3702)

file.uniprot<-"../../Blast2GO_Data/00-Downloaded/uniprot_sprot.dat.gz"
uniprot.ori<-fread(file.uniprot,fill=T,sep=" ",header = F,strip.white = F)
uniprot.ori[,Tag:=sub(" +.*","",V1)]
myMessage("The annotation file loaded.")

Tag.select<-c("ID","AC","OX","DR")
uniprot.trans<-uniprot.ori[Tag %in% Tag.select]
uniprot.trans[,Content:=substring(V1,6,length(V1))]

ID.postions<-which(uniprot.trans$Tag=="ID")
ID.lengths<-c(ID.postions[-1],nrow(uniprot.trans)+1)-ID.postions
uniprot.trans[,Num:=rep(1:length(ID.postions),ID.lengths)]
uniprot.trans<-dcast(uniprot.trans, Num ~ Tag, fun.aggregate = function(x){paste0(x,collapse = "!")},value.var = "Content")

uniprot.trans<-data.table(separate_rows(uniprot.trans,DR,sep = "!"))
# uniprot
uniprot.trans[,`:=`(DR_Tag=sub(";.*","",DR),
                    Num=Num,
                    ID=sub(";.*","",AC),
                    Accession=sub(" .*","",ID),
                    Taxid=sub(".*=([0-9]+)[; ].*","\\1",OX)
                    )]
uniprot.At<-uniprot.trans[Taxid %in% At_ids$taxid]
uniprot.At.1<-dcast(uniprot.At, ID + Accession ~ DR_Tag, fun.aggregate = function(x){paste0(x,collapse = "!")},value.var = "DR")
uniprot.At.1<-data.table(separate_rows(uniprot.At.1,GO,sep = "!"))
uniprot.At.1<-data.table(separate_rows(uniprot.At.1,TAIR,sep = "!"))
t<-head(uniprot.At,1000)
# uniprot.go<-uniprot.trans[DR_Tag=="GO",.(Num=Num,
#                                          ID=sub(";.*","",AC),
#                                          Accession=sub(" .*","",ID),
#                                          Taxid=sub(".*=([0-9]+)[; ].*","\\1",OX),
#                                          GO=DR
# )]

stop()

file.uniprot<-"../../Blast2GO_Data/01-Parsed/uniprot_sprot.tsv"
if (!file.exists(file.uniprot)) update_uniprot()
uniprot.all<-fread(file.uniprot)
uniprot.all<-fix_alternative(uniprot.all)
uniprot.At<-uniprot.all[Taxid %in% At_ids$taxid]
uniprot.At.clean<-remove_redundancy(uniprot.At,trSize = 1000,select.by =c("ID","Evidence"))
