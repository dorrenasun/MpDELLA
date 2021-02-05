#!/bin/bash
# Define a timestamp function
timestamp() {
  date +"%Y-%m-%d %H:%M:%S" >> $1 2>&1 
}



parafiles=./Parameters/*.tsv
log=blast2go.log

echo "Blast2go started." > $log
timestamp $log

for input in $parafiles; do

echo $input
Rscript Blast2GO-original.R $input >>$log 2>&1
#Rscript test.R $input >>$log 2>&1

done

echo "Blast2go ended." >> $log
timestamp $log
