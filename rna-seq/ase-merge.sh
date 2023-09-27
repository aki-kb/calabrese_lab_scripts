
######### Bash script to reformat and merge ase RNA-Seq counts from multiple datasets #########

## USAGE: sh ase-merge.sh $1
	## $1 = rootname of final outfile
	## Requirements:
		## 'samples.txt' file - see below.
		## run from top directory with a subdirectory called 'ase' with ase_analyzer11.pl outfile(s) whose rootname(s) match samples.txt

## -- ##

## create a 'samples.txt' file with file names - rootnames only, one file per line.
samplefiles=$(<samples.txt)

## create directories for outfiles
mkdir ase/sorted_chr_start
mkdir ase/sorted_chr_start_SNheader

## -- ##

## 1. sort ase_analyzer10.pl outfiles.
for file in ${samplefiles[*]}
do
	## 1a. sort ase_analyzer10.pl outfile by chr and start position (ascending).
	( head -n 1 ase/"${file}"_mm10-allelic_ase-genes ; tail -n +2 ase/"${file}"_mm10-allelic_ase-genes | sort -k2,2 -k3,3n ) > ase/sorted_chr_start/"${file}"_mm10-allelic_ase-genes_sorted
done

## 2. reformat the sorted ase-genes files.
for file in ${samplefiles[*]}
do
	## 2a. isolate gene columns
	awk 'OFS="\t" {print $1, $2, $3, $4, $5}' ase/sorted_chr_start/${file}_mm10-allelic_ase-genes_sorted > ase/starting_columns

	## 2b. isolate gene names from gene id column.
	tail -n +2 ase/starting_columns | awk -F'_' '{print $1}' > ase/gene_column

	## 2c. add header to gene names column.
	( echo "gene"; cat ase/gene_column ) > ase/gene_column2

	## 2d. isolate id, chr, start, end, strand headers.
	( printf "id\tchr\tstart\tend\tstrand\n"; tail -n +2 ase/starting_columns ) > ase/starting_columns2

	## 2e. isolate just gene id column.
	awk '{print $1}' ase/starting_columns2 > ase/id_column

	## 2f. isolate chr, start, end, strand columns.
	awk 'OFS="\t" {print $2, $3, $4, $5}' ase/starting_columns2 > ase/chr_start_end_strand_columns

	## 2g. paste together columns in the order of: id, chr, start, end, strand.
	paste ase/id_column ase/gene_column2 ase/chr_start_end_strand_columns > ase/starting_columns3

	## 2h.  create final header file
	printf "id\tchr\tstart\tend\tstrand\t${file}_B6_exon\t${file}_B6_ejc\t${file}_B6_intron\t${file}_B6_total\t${file}_cast_exon\t${file}_cast_ejc\t${file}_cast_intron\t${file}_cast_total\n" > ase/SN_header
    
	## 2i. copy headers to sorted ase-genes files.
	( cat ase/SN_header; tail -n +2 ase/sorted_chr_start/${file}_mm10-allelic_ase-genes_sorted ) > ase/sorted_chr_start_SNheader/${file}_mm10-allelic_ase-genes_sorted_SNheader
	
	## 2j. Clean this directory up, dammnit.
	rm ase/id_column
	rm ase/starting_columns
	rm ase/starting_columns2
	rm ase/chr_start_end_strand_columns
	rm ase/gene_column
	rm ase/gene_column2
	rm ase/SN_header
done

## 3. create a starting outfile base to write ase counts of all datasets.
cp ase/starting_columns3 ase/$1_mm10_ase-genes_summary.txt

## 4. add allelic counts to final outfile.
for file in ${samplefiles[*]}
do
	## 4a. isolate cast_total and b6_total read counts from ase_process.sh outfile and write to a temp.txt.
	awk 'OFS="\t" {print $13, $9}' ase/sorted_chr_start_SNheader/${file}_mm10-allelic_ase-genes_sorted_SNheader > ase/temp

	## 4b. paste together columns from final outfile and temp.txt to temp2.txt.
	paste ase/$1_mm10_ase-genes_summary.txt ase/temp > ase/temp2

	## 4c. rewrite final outfile with columns from temp2.txt.
	cat ase/temp2 > ase/$1_mm10_ase-genes_summary.txt
done

## 5. Remove unnecessary files/directories.
rm ase/temp
rm ase/temp2
rm ase/starting_columns3
rm -r ase/sorted_chr_start
rm -r ase/sorted_chr_start_SNheader

echo "YOU DID IT. woo."