
######### Pipeline to parse allelic reads and count over 10-kb bins #########

## USAGE: $ sh allelic-mm10_B6-CAST.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) create 'fastq' directory containing fastq file(s) whose rootname(s) match samples

### -- ###

## IMPORTANT (A): Before running bash script, manually create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<select.txt)

## Load modules
module load perl

## Point to genomes, lists, and scripts
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

###############################################################################

## For each sample, run the for loop below.
for file in ${sampleFiles[*]}
do
	## Parse mm10 and mm10-CAST MAPQ>=30 reads and extract only B6 and CAST SNP-overlapping reads (must overlap at least 1 SNP in the read)
	jid1=$(sbatch --time=144:00:00 --mem=12G -J ${file}.allelic -o reports/${file}.%j.allelic.out -e reports/${file}.%j.allelic.err --wrap="perl ${assets}/intersect_reads_snps18.pl sam_files/${file}_mm10_q30.sam sam_files/${file}_castmm10_q30.sam ${assets}/sanger_mm10_cast n allelic/${file}_mm10-allelic")
	jid1=$(echo $jid1 | awk '{print $NF}')

	## Count the number of allelic reads within bins across the length of eachc chromosome for tiling density analysis: 10-kb bins/no slide
	jid2=$(sbatch --dependency=afterok:$jid1 --mem=2G --time=02:00:00 -J ${file}.tiling -o reports/${file}.%j.tiling.out -e reports/${file}.%j.tiling.err --wrap="perl ${assets}/ase_analyzer8_hDbed.pl allelic/${file}_mm10-allelic_final ${assets}/all-chr_mm10_10kb-bin.bed binned_counts/${file}_mm10-allelic_10kb-bins.txt")
	jid2=$(echo $jid2 | awk '{print $NF}')
done