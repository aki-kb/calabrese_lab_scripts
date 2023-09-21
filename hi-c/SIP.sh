#!/bin/sh

#SBATCH -J sip
#SBATCH -p general
#SBATCH --mem=8G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1200
#SBATCH -o %j.sip.out
#SBATCH -e %j.sip.err

# Must be in current directory containing inter.hic and inter_30.hic files.
# Usage: sbatch SIP.sh <rootname-of-file>
        # $1 = rootname of .hic file

# ## Set the following variables to work with your system
# # export PATH=/gpfs0/biobuild/biobuilds-2016.11/bin:$PATH
#         unset MALLOC_ARENA_MAX
#         load_bwa="module load bwa"
#         load_gpu="CUDA_VISIBLE_DEVICES=0,1,2,3"
#         # Juicer directory, contains scripts/, references/, and restriction_sites/
#         # can also be set in options via -D
#         juiceDir="/proj/calabrlb/users/Keean/HiC/juicer/"
#         # default queue, can also be set in options
#         queue="general"
#         queue_time="1440"
#         # default long queue, can also be set in options
#         long_queue="general"
#         long_queue_time="10080"

# # Aiden Lab specific check
# isRice=$(hostname | awk '{if ($1~/rice/){print 1}else {print 0}}')
# isBCM=$(hostname | awk '{if ($1~/bcm/){print 1}else {print 0}}')
# isVoltron=0
# ## path additionals, make sure paths are correct for your system
# ## use cluster load commands
# if [ $isRice -eq 1 ]
# then
#     	myPath=/bin:$PATH
#         load_bwa="module load BioBuilds/2015.04"
#         load_java="module load Java/8.0.3.22"
#         load_gpu="module load gcccuda/2016a;module load CUDA/8.0.54;"
#         # Juicer directory, contains scripts/, references/, and restriction_sites/
#         # can also be set in options via -D
#         juiceDir="/projects/ea14/juicer" ### RICE
#         # default queue, can also be set in options via -q
#         queue="general"
#         queue_time="1440"
#         # default long queue, can also be set in options via -l
#         long_queue="general"
#         long_queue_time="1440"
# elif [ $isBCM -eq 1 ]
# then
#     	# Juicer directory, contains scripts/, references/, and restriction_sites/
#         # can also be set in options via -D
#         juiceDir="/storage/aiden/juicer/"
#         # default queue, can also be set in options via -q
#         queue="general"
#         queue_time="1200"
#         # default long queue, can also be set in options via -l
#         long_queue="general"
#         long_queue_time="3600"
# else
#     	isVoltron=1
#         # export PATH=/gpfs0/biobuild/biobuilds-2016.11/bin:$PATH
#         unset MALLOC_ARENA_MAX
#         load_bwa="module load bwa"
#         load_gpu="CUDA_VISIBLE_DEVICES=0,1,2,3"
#         # Juicer directory, contains scripts/, references/, and restriction_sites/
#         # can also be set in options via -D
#         juiceDir="/proj/calabrlb/users/Keean/HiC/juicer/"
#         # default queue, can also be set in options
#         queue="general"
#         queue_time="1440"
#         # default long queue, can also be set in options
#         long_queue="general"
#         long_queue_time="10080"
# fi

# ${load_java}

module load java

HiC="/proj/calabrlb/users/Keean/HiC"

## make output directory
#mkdir $1_B6_chr17_SIP
#mkdir $1_CASTEiJ_chr17_SIP

## run SIP on .hic file; default is 5kb (recommended for high resolution data)
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic ./inter.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes $1_SIP ${HiC}/juicer/scripts/juicer_tools.jar
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic ./inter_30.hic ${HiC}/Mus_musculus_assembly9_norandom.chrom.sizes $1_SIP ${HiC}/juicer/scripts/juicer_tools.jar

## run SIP on .hic files; loops calls for 5, 10, and 25 kb resolutions
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic ./$1.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes $1_SIP ${HiC}/juicer/scripts/juicer_tools.jar -res 5000 -factor 4
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic ./$1_q30.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes $1_q30_SIP ${HiC}/juicer/scripts/juicer_tools.jar -res 5000 -factor 4

## as of 15.04.2021 - parameters provided by Eric (tested by Katie Reed)
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes $1_SIP.out ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_q30.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes $1_q30_SIP.out ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05

# Mus_musculus_assembly9_norandom.chrom.sizes2 excludes chrY and chrM, which allelic data does not hold - it will error out otherwise.
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_B6_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes2 $1_B6_chr17_SIP ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_CASTEiJ_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom.sizes2 $1_CASTEiJ_chr17_SIP ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05

# Mus_musculus_assembly9_norandom.chrom17.sizes.
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_B6_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom17.sizes $1_B6_chr17_SIP ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05
#java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_CASTEiJ_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom17.sizes $1_CASTEiJ_chr17_SIP ${HiC}/juicer/scripts/juicer_tools.jar -factor 4 -g 2.0 -t 2000 -fdr 0.05

java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_B6_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom17.sizes $1_B6_chr17_SIP_50kb ${HiC}/juicer/scripts/juicer_tools.jar -res 50000 -g 2.0 -t 2000 -fdr 0.05
java -jar ${HiC}/SIP/SIP_HiC_v1.5.jar hic $1_CASTEiJ_chr17.hic ${HiC}/juicer/references/Mus_musculus_assembly9_norandom.chrom17.sizes $1_CASTEiJ_chr17_SIP_50kb ${HiC}/juicer/scripts/juicer_tools.jar -res 50000 -g 2.0 -t 2000 -fdr 0.05

