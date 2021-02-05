#!/bin/bash
# head -n 10 Blast_Mpa_v3.1-uniprot.tsv > test.tsv
# Check how many cpu cores
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	ncpu=$(nproc --all)
elif [[ "$OSTYPE" == "darwin"* ]]; then
	ncpu=$(getconf _NPROCESSORS_ONLN)
fi

fasta_query=query.fasta
fasta_uniprot=uniprot_selected.fasta

dir_results=Blast1_results

logfile=${dir_results}/Blast1_log.txt
# blast_out0=${dir_results}/Blast_query-self.tsv
blast_out1=${dir_results}/Blast_query-uniprot.tsv


mkdir $dir_results
echo "Number of cores: $ncpu" >> $logfile 2>&1
now=$(date +"%Y-%m-%d %T")
echo "\nBlast1 started: $now" >> $logfile 2>&1

# echo "\nBlast round one: $fasta_query self" >> $logfile 2>&1
# makeblastdb -in $fasta_query -parse_seqids -dbtype prot -out $fasta_query >> $logfile 2>&1
# blastp -db $fasta_query -query $fasta_query -out $blast_out0 -evalue 0.001 -outfmt "6 std ppos" -num_threads $ncpu
# rm ${fasta_query}.*

echo "\nBlast round two: $fasta_query vs $fasta_uniprot" >> $logfile 2>&1
makeblastdb -in $fasta_uniprot -parse_seqids -dbtype prot -out $fasta_uniprot >> $logfile 2>&1
blastp -db $fasta_uniprot -query $fasta_query -out $blast_out1 -evalue 0.001 -outfmt "6 std ppos" -num_threads $ncpu
rm ${fasta_uniprot}.*


now=$(date +"%Y-%m-%d %T")
echo "\nBlast1 finished: $now" >> $logfile 2>&1

tar -zcvf ${dir_results}.tar.gz $dir_results

gsutil cp ${dir_results}.tar.gz gs://blast2go/
sudo poweroff

#Notes: how to run on Google cloud
# gsutil cp gs://blast2go/Blast1.zip ./
# unzip Blast1.zip -d Blast1
# cd Blast1
# setsid ./Blast1.sh &
# top ## to check if blast is running