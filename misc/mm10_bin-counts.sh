
######### Total, non-allelic ChIP-Seq pipeline for counting mm10 reads in 10-kb bins #########

## USAGE: $ sh mm10_bin-counts.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) create 'fastq' directory containing fastq file(s) whose rootname(s) match samples

### -- ###

## IMPORTANT (A): Before running bash script, manually create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<samples.txt)

## Load modules
module load bedtools

## Point to genomes, lists, and scripts
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

###############################################################################

## For each sample, run the for loop below.
for file in ${sampleFiles[*]}
do

	## Count the number of mm19 reads in genomic bins for tiling density analysis: 10-kb bins/no slide
	sbatch --mem=2G --time=04:00:00 -J ${file}.tiling -o reports/${file}.%j.tiling.out -e reports/${file}.%j.tiling.err --wrap="bedtools coverage -counts -sorted -g ${assets}/mm10.chr.sizes.sorted -a ${assets}/all-chr_mm10_10kb-bin_sorted.bed -b bam_files/${file}_mm10_q30.sorted.bam > binned_counts/${file}_mm10_all-chr_10kb-bins.txt"

done