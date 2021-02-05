##Downloaded
file.GO<-"../../Blast2GO_Data/00-Downloaded/go.obo"
# file.fasta.query<-"../../Blast2GO_Data/00-Downloaded/Mpolymorphav3.1.allTrs.pep_annot.fa"
file.fasta.ver5<-"../../Blast2GO_Data/00-Downloaded/MpTak1v5.1_r1.protein.fasta"
file.fasta.ver3<-"../../Blast2GO_Data/00-Downloaded/Mpolymorphav3.1.allTrs.pep_annot.fa.gz"



file.sprot.dat<-"../../Blast2GO_Data/00-Downloaded/uniprot_sprot.dat.gz"
file.sprot.fasta<-"../../Blast2GO_Data/00-Downloaded/uniprot_sprot.fasta.gz"
file.sprot.reldate<-"../../Blast2GO_Data/00-Downloaded/uniprot_sprot_date.txt"

file.noiea.gaf<-"../../Blast2GO_Data/00-Downloaded/goa_uniprot_all_noiea.gaf.gz"
file.noiea.reldate<-"../../Blast2GO_Data/00-Downloaded/goa_uniprot_all_noiea_summary.txt"

file.taxdmp<-"../../Blast2GO_Data/00-Downloaded/taxdmp.zip"

##Parsed
file.query.parsed<-"../../Blast2GO_Data/01-Parsed/Query/Mp_genome_v3.1.tsv"
file.fasta.ver5f<-"../../Blast2GO_Data/01-Parsed/Query/Mp_genome_v5+female.fasta"

file.taxon<-"../../Blast2GO_Data/01-Parsed/Taxonmy/myTaxonomy.tsv"
file.unknown<-"../../Blast2GO_Data/01-Parsed/Taxonmy/Manual_taxids.tsv"
file.plantid<-"../../Blast2GO_Data/01-Parsed/Taxonmy/Plant_IDs.tsv"
file.nonplant.BP<-"../../Blast2GO_Data/01-Parsed/Taxonmy/NonPlants_GOs_BP.tsv"

file.uniprot.go<-"../../Blast2GO_Data/01-Parsed/Uniprot/uniprot_go.tsv"
file.uniprot.go.selected<-"../../Blast2GO_Data/01-Parsed/Uniprot/uniprot_go_filtered.tsv"

file.uniprot.download<-"../../Blast2GO_Data/01-Parsed/Uniprot/downloaded.gb"
file.uniprot.fasta<-"../../Blast2GO_Data/01-Parsed/Uniprot/uniprot_selected.fasta"
file.uniprot.nofasta<-"../../Blast2GO_Data/01-Parsed/Uniprot/uniprot_no_info.tsv"

#Blast1
path.Blast1<-"../../Blast2GO_Data/01-Parsed/Blast1/"
file.Blast1.query<-paste0(path.Blast1,"query.fasta")
file.Blast1.zip<-paste0(path.Blast1,"Blast1.zip")
file.Blast1.result.self<-paste0(path.Blast1,"Blast1_results/Blast_query-self.tsv")
file.Blast1.result.uniprot<-paste0(path.Blast1,"Blast1_results/Blast_query-uniprot.tsv")

#Blast2
path.Blast2<-"../../Blast2GO_Data/01-Parsed/Blast2/"
file.Blast2.query<-paste0(path.Blast2,"query.fasta")
file.Blast2.hits<-paste0(path.Blast2,"uniprot_hits.fasta")
file.Blast2.zip<-paste0(path.Blast2,"Blast2.zip")
file.Blast2.groups<-paste0(path.Blast2,"Groups.tsv")
path.Blast2.results<-paste0(path.Blast2,"Blast2_results/")
file.Blast2.result.query<-paste0(path.Blast2.results,"Blast_hits-query.tsv")

#OrthoFinder
path.ortho<-"../../Blast2GO_Data/01-Parsed/OrthoFinder/"

#Blast2GO
path.blast2go<-"../../Blast2GO_Data/02-Results/Blast2GO/"

