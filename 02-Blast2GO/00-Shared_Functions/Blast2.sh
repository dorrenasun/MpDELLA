#!/bin/bash
# head -n 10 Blast_Mpa_v3.1-uniprot.tsv > test.tsv
# Check how many cpu cores
fasta_query=query.fasta
fasta_hits=uniprot_hits.fasta

dir_results=./Blast2_results
mkdir $dir_results
blast_out2=${dir_results}/Blast_hits-query.tsv

logfile=${dir_results}/log.txt


if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	ncpu=$(nproc --all)
elif [[ "$OSTYPE" == "darwin"* ]]; then
	ncpu=$(getconf _NPROCESSORS_ONLN)
fi
echo "Number of cores: $ncpu" >> $logfile 2>&1

now=$(date +"%Y-%m-%d %T")
echo "Blast2 started: $now" >> $logfile 2>&1
echo "\nBlast round three: $fasta_hits vs $fasta_query" >> $logfile 2>&1
makeblastdb -in $fasta_query -parse_seqids -dbtype prot -out $fasta_query >> $logfile 2>&1
blastp -db $fasta_query -query $fasta_hits -out $blast_out2 -evalue 0.001 -outfmt "6 std ppos" -num_threads $ncpu
rm ${fasta_query}.*


fastalist=Group*.fasta

for f in $fastalist
do
	echo "Processing $f" >> $logfile 2>&1
	makeblastdb -in $f -parse_seqids -dbtype prot -out $f >> $logfile 2>&1
	outfile=${dir_results}/Blast-$f.tsv
	# echo $outfile
	blastp -db $f -query $f -out $outfile -evalue 0.001 -outfmt "6 std ppos" -num_threads $ncpu
	rm $f.*
done

now=$(date +"%Y-%m-%d %T")
echo "Blast2 finished: $now" >> $logfile 2>&1


tar -zcvf ${dir_results}.tar.gz $dir_results

gsutil cp ${dir_results}.gz gs://blast2go/
sudo poweroff
