
######### Filter mm10-aligned reads by MAPQ>=30. #########

## USAGE: $ sh sam_q30-filter.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) sam_files/*_mm10.sam files must exist


## IMPORTANT (A): Before running bash script, manually create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<samples.txt)

## Load modules
module load samtools

## Point to genomes, lists, and scripts
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

###############################################################################

## For each sample, run the for loop below.
for file in ${sampleFiles[*]}
do

	## Filter mm10- and mm10-CAST reads for only MAPQ>=30
	sbatch --time=04:00:00 --wrap="samtools view -q 30 sam_files/${file}_mm10.sam > sam_files/${file}_mm10_q30.sam"

done