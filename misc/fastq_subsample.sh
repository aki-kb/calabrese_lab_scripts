## IMPORTANT (A): Before running bash script, manually create a 'select.txt' file with file names - rootnames only, one file per line.
sampleFiles=$(<select.txt)

## Load modules
module load seqtk

for file in ${sampleFiles[*]}
do
    if [ -f fastq/${file}.fastq.gz ]
	then
		echo "${file}.fastq.gz exists, gunzip-ing now."
		gunzip fastq/${file}.fastq.gz
		echo "Gunzip complete, proceeding with subsampling."
	else
		echo "${file}.fastq exists, proceeding with subsampling."
	fi

    # Randomly subsample 25 million reads from FASTQ.
    sbatch -p general --time=04:00:00 --mem=24G -J ${file}.subsample -o reports/${file}.%j.subsample.out -e reports/${file}.%j.subsample.err --wrap="seqtk sample -s 123 fastq/${file}.fastq 25000000 > fastq/${file}_25M.fastq"

done