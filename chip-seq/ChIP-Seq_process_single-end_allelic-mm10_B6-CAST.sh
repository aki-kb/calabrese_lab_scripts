
######### Allelic ChIP-Seq pipeline for single-end reads against mm10 and CAST-mm10 genome builds #########

## USAGE: $ sh ChIP-seq_process_single-end_allelic-mm10_B6-CAST.sh
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
	wiggles
	allelic
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

	## 2) Align fastq against mm10 (B6) and mm10-CAST genomes with bowtie2
		## Make sure to use bowtie2 indexed genomes associated with bowtie2 version 2.4.5 and genomes
	jid1=$(sbatch -p general -N 1 -n 24 --time=4:00:00 --mem=10G -J ${file}.mm10.bowtie2 -o reports/${file}.%j.mm10.bowtie2.out -e reports/${file}.%j.mm10.bowtie2.err --wrap="bowtie2 -q -p 24 -x ${assets}/Bowtie2.4.5_genome_index_mm10/genome -U fastq/${file}.fastq -S sam_files/${file}_mm10.sam")
	jid1=$(echo $jid1 | awk '{print $NF}')
	jid2=$(sbatch -p general -N 1 -n 24 --time=4:00:00 --mem=10G -J ${file}.castmm10.bowtie2 -o reports/${file}.%j.castmm10.bowtie2.out -e reports/${file}.%j.castmm10.bowtie2.err --wrap="bowtie2 -q -p 24 -x ${assets}/Bowtie2.4.5_genome_index_CASTmm10/genome -U fastq/${file}.fastq -S sam_files/${file}_castmm10.sam")
  	jid2=$(echo $jid2 | awk '{print $NF}')

	## 3) Filter mm10- and mm10-CAST reads for only MAPQ>=30
	jid3=$(sbatch --dependency=afterok:$jid1 --time=04:00:00 --wrap="samtools view -q 30 sam_files/${file}_mm10.sam > sam_files/${file}_mm10_q30.sam")
	jid3=$(echo $jid3 | awk '{print $NF}')
	jid4=$(sbatch --dependency=afterok:$jid2 --time=04:00:00 --wrap="samtools view -q 30 sam_files/${file}_castmm10.sam > sam_files/${file}_castmm10_q30.sam")
	jid4=$(echo $jid4 | awk '{print $NF}')

	## 4) Parse mm10 and mm10-CAST MAPQ>=30 reads and extract only B6 and CAST SNP-overlapping reads (must overlap at least 1 SNP in the read)
	jid5=$(sbatch --dependency=afterok:$jid3,$jid4 --time=72:00:00 --mem=10G -J ${file}.allelic -o reports/${file}.%j.allelic.out -e reports/${file}.%j.allelic.err --wrap="perl ${assets}/intersect_reads_snps18.pl sam_files/${file}_mm10_q30.sam sam_files/${file}_castmm10_q30.sam ${assets}/sanger_mm10_cast n allelic/${file}_mm10-allelic")
	jid5=$(echo $jid5 | awk '{print $NF}')

	## 5) Generate a wiggle file from the mm10 MAPQ>=30 reads for UCSC genome browser viewing
	jid6=$(sbatch --dependency=afterok:$jid3 --time=96:00:00 --mem=2G -J ${file}.wiggle -o reports/${file}.%j.wiggle.out -e reports/${file}.%j.wiggle.err --wrap="perl ${assets}/bigsam_to_wig_mm10_wcigar4.pl sam_files/${file}_mm10_q30.sam ${assets}/mm10.chr.sizes wiggles/${file}_mm10 purple y n 50")
	jid6=$(echo $jid6 | awk '{print $NF}')

	## 6) Count the number of allelic reads within bins across the length of eachc chromosome for tiling density analysis: 10-kb bins/no slide
	jid7=$(sbatch --dependency=afterok:$jid5 --mem=2G --time=02:00:00 -J ${file}.tiling -o reports/${file}.%j.tiling.out -e reports/${file}.%j.tiling.err --wrap="perl ${assets}/ase_analyzer8_hDbed.pl allelic/${file}_mm10-allelic_final ${assets}/all-chr_mm10_10kb-bin.bed binned_counts/${file}_mm10-allelic_10kb-bins.txt")
	jid7=$(echo $jid7 | awk '{print $NF}')

	## 7) Gzip fastq
	jid8=$(sbatch --dependency=afterok:$jid1,$jid2 --wrap="gzip fastq/${file}.fastq")
	jid8=$(echo $jid8 | awk '{print $NF}')
done

## Notes on dependencies:
	## jid1 and jid2 are independent of each other
	## jid3 and jid4 are independent of each other, but dependent on jid1 and jid2, respectively
	## jid5 is dependent on jid3 and jid4
	## jid6 is dependent on jid3
	## jid7 is dependent on jid5
	## jid8 is dependent on jid1 and jid2

## Notes on output files:
	## sam_files/ - contains all sam files generated from bowtie2 alignments
	## wiggles/ - contains all wiggle files generated from sam files
		## wiggles can be scaled scaled separately or gzip-ed for UCSC genome browser viewing
	## allelic/ - contains all allelic files generated from sam files
	## binned_counts/ - contains all binned counts files generated from allelic files
		## binned counts files can be used for tiling density analysis in RStudio - see tiling_density_analysis.Rmd
	## reports/ - contains all .out and .err files generated from sbatch jobs
		## *bowtie2.out files contain alignment statistics, including the number of total reads for RPM calculations