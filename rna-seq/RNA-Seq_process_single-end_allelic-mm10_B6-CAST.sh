######### Allelic RNA-Seq pipeline for single-end reads against mm10 and CAST-mm10 genome builds #########

## USAGE: sh RNA-Seq_process_single-end_allelic-mm10_B6-CAST.sh
        ## Input REQUIREMENTS:
			## A) create 'samples.txt' file - list of rootnames for files (separated by \n) to be processed
			## B) create 'fastq' directory containing fastq file(s) whose rootname(s) match samples

### -- ###

## load modules
module load star/2.7.9a
module load samtools
module load perl

## create a 'samples.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<samples.txt)

## create directories for output files
mkdir sam_files
mkdir wiggles
mkdir allelic
mkdir ase
mkdir reports

## create alias for directory pathways
assets="/proj/calabrlb/users/Keean/bioinformatics/calabrese_lab_scripts/assets"

## for each file as specified in 'samples.txt', apply below.
for file in ${sampleFiles[*]}
do
	## 1. Check status of fastq; fastq needs to not be gzipped for alignments and downstream processing.
	if [[ -f fastq/${file}.fastq.gz ]]
	then
		echo "${file}.fastq.gz exists, gunzip-ing now."
		gunzip fastq/${file}.fastq.gz
		echo "Gunzip complete, proceeding with alignments."
	else
		echo "${file}.fastq exists, proceeding with alignments."
	fi

	## 2. Map fastq against mm10 B6 and CAST genomes with STAR; SINGLE-END MODE
	# Make sure to use STAR indexed genomes associated with specific STAR versions and genomes - directories must be in working directory
	# STAR alignments only considered uniquely mapping reads (q=255); "--outFilterMultimapNmax 1" passes read to alignment outfile only if aligned to 1 locus
	jid1=$(sbatch -p general --time=02:00:00 -N 1 -n 12 --mem=32g -J ${file}.star.mm10 -o reports/${file}.%j.star.mm10.out -e reports/${file}.%j.star.mm10.err --wrap="STAR --genomeDir ${assets}/STARv2.7.9a_genome_index_mm10 --runThreadN 12 --outFilterMultimapNmax 1 --readFilesIn fastq/${file}.fastq --outFileNamePrefix ${file}_mm10_ --outSAMtype SAM")
	jid1=$(echo $jid1 | awk '{print $NF}')
	jid2=$(sbatch -p general --time=02:00:00 -N 1 -n 12 --mem=32g -J ${file}.star.castmm10 -o reports/${file}.%j.star.castmm10.out -e reports/${file}.%j.star.castmm10.err --wrap="STAR --genomeDir ${assets}/STARv2.7.9a_genome_index_CASTmm10 --runThreadN 12 --outFilterMultimapNmax 1 --readFilesIn fastq/${file}.fastq --outFileNamePrefix ${file}_castmm10_ --outSAMtype SAM")
  	jid2=$(echo $jid2 | awk '{print $NF}')

	## 3. Separate reverse-stranded mm10-aligned data by strand for wiggle generator
	jid3=$(sbatch --dependency=afterok:$jid1 --time=04:00:00 --wrap="samtools view -F 0x10 ${file}_mm10_Aligned.out.sam > sam_files/${file}_mm10_neg.sam")
	jid3=$(echo $jid3 | awk '{print $NF}')
	jid4=$(sbatch --dependency=afterok:$jid1 --time=04:00:00 --wrap="samtools view -f 0x10 ${file}_mm10_Aligned.out.sam > sam_files/${file}_mm10_pos.sam")
	jid4=$(echo $jid4 | awk '{print $NF}')

	## 4. Prepare wiggle files for mm10/B6-aligned data
	jid5=$(sbatch --dependency=afterok:$jid3 --mem=2G --time=16:00:00 -J ${file}.wiggle.neg -o reports/${file}.%j.wiggle.neg.out -e reports/${file}.%j.wiggle.neg.out --wrap="perl ${assets}/bigsam_to_wig_mm10_wcigar4.pl sam_files/${file}_mm10_neg.sam ${assets}/mm10.chr.sizes wiggles/${file}_mm10_neg orange y n 50")
	jid5=$(echo $jid5 | awk '{print $NF}')
	jid6=$(sbatch --dependency=afterok:$jid4 --mem=2G --time=16:00:00 -J ${file}.wiggle.pos -o reports/${file}.%j.wiggle.pos.out -e reports/${file}.%j.wiggle.pos.out --wrap="perl ${assets}/bigsam_to_wig_mm10_wcigar4.pl sam_files/${file}_mm10_pos.sam ${assets}/mm10.chr.sizes wiggles/${file}_mm10_pos black y n 50")
	jid6=$(echo $jid6 | awk '{print $NF}')

	## 5. Parse and distributing B6 and CAST single-end reads; SINGLE-END MODE
	jid7=$(sbatch --dependency=afterok:$jid1,$jid2 -p general --time=48:00:00 --mem=32g -J ${file}.intersect -o reports/${file}.%j.intersect.out -e reports/${file}.%j.intersect.err --wrap="perl ${assets}/intersect_reads_snps18.pl ${file}_mm10_Aligned.out.sam ${file}_castmm10_Aligned.out.sam ${assets}/sanger_mm10_cast n allelic/${file}_mm10-allelic")
	jid7=$(echo $jid7 | awk '{print $NF}')

	## 6. Gzip wiggle files, upload wig.gz files to UCSC
	jid8=$(sbatch --dependency=afterok:$jid1,$jid2 -J ${file}.gzip.fastq -o reports/${file}.%j.gzip.fastq.out -e reports/${file}.%j.gzip.fastq.err --wrap="gzip fastq/${file}.fastq")
	jid8=$(echo $jid8 | awk '{print $NF}')

	## 7. Count allelic reads over mm10 gene features (gencode.vM25.annotation.gtf): SINGLE-END MODE
	jid9=$(sbatch --dependency=afterok:$jid7 -p general --time=12:00:00 --mem=16g -J ${file}.ase11 -o reports/${file}.%j.ase11.out -e reports/${file}.%j.ase11.err --wrap="perl ${assets}/ase_analyzer11.pl ${assets}/gencode.vM25.annotation.gtf allelic/${file}_mm10-allelic_final n ase/${file}_mm10-allelic_ase-genes")
	jid9=$(echo $jid9 | awk '{print $NF}')

	## 8. Move sam files to sam_files directory
	jid10=$(sbatch --dependency=afterok:$jid5 -p general -J ${file}.mvsam --wrap="mv ${file}_mm10_Aligned.out.sam ${file}_castmm10_Aligned.out.sam sam_files/")
	jid10=$(echo $jid10 | awk '{print $NF}')
done