# In house Blast2GO annotation of the Marchantia ver5.1 genome

In this analysis, functional annotations were generated for the version 5.1 genome of Marchantia polymorpha using an _in house_ Blast2GO analysis.

## What is Bast2GO
Please see [the official site] (https://www.blast2go.com) for the original software, and 
> Götz S., Garcia-Gomez JM., Terol J., Williams TD., Nagaraj SH., Nueda MJ., Robles M., Talon M., Dopazo J. and Conesa A. (2008). High-throughput functional annotation and data mining with the Blast2GO suite. Nucleic acids research, 36(10), 3420-35. doi:https://doi.org/10.1093/nar/gkn176 
for the description of the original algorithm.

## How to run the script
Download this folder with all subdirectories. Run "Blast2GO-original.R" from the "01-Blast2GO-R" directory.

## Target sequences
All transcripts from the male MpTak1v5.1 genome, plus female-specific transcripts from v3.1 (scaffolds 0017, 0018, 0210, 0227, 0240, 0250 and 0497) were retrieved from https://marchantia.info/ and used in this analysis.

## Blast database
Fasta sequences were retrieved from the IDs included below and used for local blastp search:
* GO annotations for manually reviewed Swiss-Prot entries, from UniProt database 
	https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/
	(Clean annotation with no qualifier flags)
* Non-IEA annotations from all entries in UniProtKB, provided by the Gene Ontology knowledgebase
	http://current.geneontology.org/annotations/goa_uniprot_all_noiea.gaf.gz 
	(GO terms with qualifier flags like 'NOT', 'contributes_to' or 'colocalizes_with' were filtered out from this dataset.)

GO annotations from both sources were combined as the mapping database.

## Local blast search
After downloading and retrieving database protein sequences, they were packed together with the target sequences and a shell script to run local blastp on any computer. This step is time- and resource-requiring, and it would be better done on a linux/unix computer or server with high performance.

## Mapping & Scoring
After blast search, top blast hits were extracted from the result, and mapped to their annotations. Scoring was done according to the original algorithm, with:

	𝑆𝑐𝑜𝑟𝑒= max(𝑆𝑖𝑚𝑖𝑙𝑎𝑟𝑖𝑡𝑦 × 𝐸𝐶𝑤) + #𝑆𝑢𝑏𝑡𝑒𝑟𝑚 × 𝐺𝑂𝑤
	* #𝑆𝑢𝑏𝑡𝑒𝑟𝑚 is the number of their offspring terms with blast hits.
	* Ecw is a weight given to annotations based on their evidence codes (reliability).

The scores were then filtered with a manually selected threshold. Different thresholds could be set for different categories of GO terms (BP, CC or MF). All parameters are defined by the file in "03-Parameters". If you want to test different combination of parameters, use the bash shell script "Run_blast2go.sh" in "01-Blast2GO-R" to do the job.

## Redundancy removal of the results
In the end, repetitive and redundant terms which could be inferred from a more specific term based on the GO hierarchy was removed.


