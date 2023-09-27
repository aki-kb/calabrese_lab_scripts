
######### Bash script to generate mm10 wiggles #########

## USAGE: $ sh mm10_wiggle.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) reports/ and sam_files/*mm10_q30.sam outputs from ChIP-seq_process_single-end_allelic-mm10_B6-CAST.sh

### -- ###

## IMPORTANT (A): Before running bash script, manually create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<samples.txt)

## Load modules
module load perl

## Point to genomes, lists, and scripts
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

###############################################################################

## For each sample, run the for loop below.
for file in ${sampleFiles[*]}
do

	## NEW script: Generate a wiggle file from the mm10 MAPQ>=30 reads for UCSC genome browser viewing
	# sbatch --time=96:00:00 --mem=2G -J ${file}.wiggle -o reports/${file}.%j.wiggle.out -e reports/${file}.%j.wiggle.err --wrap="perl ${assets}/bigsam_to_wig_mm10_wcigar4.pl sam_files/${file}_mm10_q30.sam ${assets}/mm10.chr.sizes wiggles/${file}_mm10 purple y n 50"

	## OLD script: Generate a wiggle file from the mm10 MAPQ>=30 reads for UCSC genome browser viewing
	sbatch --time=12:00:00 --mem=2G -J ${file}.wiggle -o reports/${file}.%j.wiggle.out -e reports/${file}.%j.wiggle.err --wrap="perl ${assets}/bigbowtie_to_wig3_mm10.pl sam_files/${file}_mm10_q30.sam wiggles/${file}_mm10 purple"

done

