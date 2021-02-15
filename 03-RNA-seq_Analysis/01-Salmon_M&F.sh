#!/bin/bash

timestamp() {
  date +"%Y-%m-%d %H:%M:%S" >> $1 2>&1 
}


topdir=../FASTQ

files_Tak=`(find "$topdir" -type f \( -name '*.fastq.gz' \) | grep Tak1)`
files_pif=`(find "$topdir" -type f \( -name '*.fastq.gz' \) | grep pif)`

rm -R Salmon
mkdir Salmon
log=./Salmon/Salmon.log
touch $log

echo "Salmon analysis started." > $log
echo "Salmon version:" >>$log 2>&1
salmon -v >> $log 2>&1
echo "--------------------------" >>$log 2>&1

timestamp $log
for input in $files_Tak; do
	output=$(awk -F/ '{ print $4 "_" $5 }' <<<"${input}")
	output=./Salmon/$(basename $output .fastq.gz)

	echo "Processing: ${input}" >> $log
	salmon quant -i ../../Mpver5.1/salmon_index_Tak1/ -l A -r $input  -p 4 --validateMappings -o $output >> $log 2>&1

done
for input in $files_pif; do
	output=$(awk -F/ '{ print $4 "_" $5 }' <<<"${input}")
	output=./Salmon/$(basename $output .fastq.gz)

	echo "Processing: ${input}" >> $log
	salmon quant -i ../../Mpver5.1/salmon_index_Female/ -l A -r $input  -p 4 --validateMappings -o $output >> $log 2>&1

done

timestamp $log