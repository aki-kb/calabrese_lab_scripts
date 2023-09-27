
######### Total, non-allelic ChIP-Seq pipeline for single-end reads against mm10 genome build #########

## USAGE: $ sh ChIP-seq_process_single-end_mm10.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) create 'fastq' directory containing fastq file(s) whose rootname(s) match samples

### -- ###

## IMPORTANT (A): Before running bash script, manually create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<samples.txt)

## Load modules
module load gcc/11.2.0
module load bowtie2/2.4.5
module load samtools
module load perl
module load bedtools

## Check to see if fastq/ directory exists; if not, aborts pipeline.
if [ -d fastq ]
then
	echo "fastq directory exists, proceeding."
else
	echo "fastq directory does not exist, aborting pipeline."
	exit 1
fi

## Create new directories for output files if they don't already exist.
newDir=(
	sam_files
	bam_files
	wiggles
	binned_counts
	reports
	)

for dir in ${!newDir[*]}
do
	if [ -d ${newDir[$dir]} ]
	then
		echo "${newDir[$dir]}/ directory exists."
	else
		mkdir ${newDir[$dir]}
		echo "${newDir[$dir]}/ directory created."
	fi
done

## Point to genomes, lists, and scripts
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

###############################################################################

## For each sample, run the for loop below.
for file in ${sampleFiles[*]}
do
	## 1) Check status of fastq; fastq needs to not be gzipped for alignments and downstream processing.
	if [ -f fastq/${file}.fastq.gz ]
	then
		echo "${file}.fastq.gz exists, gunzip-ing now."
		gunzip fastq/${file}.fastq.gz
		echo "Gunzip complete, proceeding with alignments."
	else
		echo "${file}.fastq exists, proceeding with alignments."
	fi

	## 2) Align fastq against mm10 genomes with bowtie2
		## Make sure to use bowtie2 indexed genomes associated with bowtie2 version 2.4.5 and genomes
	jid1=$(sbatch -p general -N 1 -n 24 --time=4:00:00 --mem=10G -J ${file}.mm10.bowtie2 -o reports/${file}.%j.mm10.bowtie2.out -e reports/${file}.%j.mm10.bowtie2.err --wrap="bowtie2 -q -p 24 -x ${assets}/Bowtie2.4.5_genome_index_mm10/genome -U fastq/${file}.fastq -S sam_files/${file}_mm10.sam")
	jid1=$(echo $jid1 | awk '{print $NF}')

	## 3) Filter and convert mm10-aligned reads by MAPQ>=30 in sam format
	jid2=$(sbatch --dependency=afterok:$jid1 --time=04:00:00 --wrap="samtools view -q 30 sam_files/${file}_mm10.sam > sam_files/${file}_mm10_q30.sam")
	jid2=$(echo $jid2 | awk '{print $NF}')

	## 4) Generate a wiggle file from the mm10 MAPQ>=30 reads for UCSC genome browser viewing
	jid3=$(sbatch --dependency=afterok:$jid2 --time=96:00:00 --mem=2G -J ${file}.wiggle -o reports/${file}.%j.wiggle.out -e reports/${file}.%j.wiggle.err --wrap="perl ${assets}/bigsam_to_wig_mm10_wcigar4.pl sam_files/${file}_mm10_q30.sam ${assets}/mm10.chr.sizes wiggles/${file}_mm10 purple y n 50")
	jid3=$(echo $jid3 | awk '{print $NF}')

	## 5) Filter and convert mm10-aligned reads by MAPQ>=30 in bam format
	jid4=$(sbatch --dependency=afterok:$jid1 -J ${file}.map30.bam --time=02:00:00 --wrap="samtools view -bhq 30 sam_files/${file}_mm10.sam > bam_files/${file}_mm10_q30.bam")
	jid4=$(echo $jid4 | awk '{print $NF}')

	## 6) Sort bam file by coordinate
	jid5=$(sbatch --dependency=afterok:$jid4 -J ${file}.bam.sorting --mem=2G --time=04:00:00 --wrap "samtools sort -o bam_files/${file}_mm10_q30.sorted.bam bam_files/${file}_mm10_q30.bam")
	jid5=$(echo $jid5 | awk '{print $NF}') 

	## 7) Count the number of mm10 reads in genomic bins for tiling density analysis: 10-kb bins/no slide
	jid6=$(sbatch --dependency=afterok:$jid5 --mem=2G --time=04:00:00 -J ${file}.tiling -o reports/${file}.%j.tiling.out -e reports/${file}.%j.tiling.err --wrap="bedtools coverage -counts -sorted -g ${assets}/mm10.chr.sizes.sorted -a ${assets}/all-chr_mm10_10kb-bin_sorted.bed -b bam_files/${file}_mm10_q30.sorted.bam > binned_counts/${file}_mm10_10kb-bins.txt")
	jid6=$(echo $jid6 | awk '{print $NF}')

	## 8) Gzip fastq
	jid7=$(sbatch --dependency=afterok:$jid1 --wrap="gzip fastq/${file}.fastq")
	jid7=$(echo $jid7 | awk '{print $NF}')

done