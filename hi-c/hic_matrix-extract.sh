#!/bin/sh

#SBATCH -J hicextract
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o %j.hicextract.out
#SBATCH -e %j.hicextract.err

# Script to extract Hi-C viewpoint counts from .hic file.

# Must be in current directory containing .hic file.
# USAGE: $ sbatch /proj/calabrlb/users/Keean/HiC/Scripts/hic_matrix-extract.sh $1 $2
        # $1 = rootname of hic file.
	# $2 = resolution (e.g. 10000)

# # Set the following variables to work with your system
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

directory="/proj/calabrlb/users/Keean/HiC/juicer/scripts"

# Comment out the viewpoints you don't want to extract.

# Airn vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:12900000:13050000 17 BP $2 $1_airn-vs-chr17_$2-res_obs_no-norm.txt
#java -jar ${directory}/juicer_tools.jar dump oe NONE $1.hic 17:12900000:13050000 17 BP $2 $1_airn-vs-chr17_$2-res_oe_no-norm.txt

# Slc22a3 vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:12612840:12700570 17 BP $2 $1_slc-gene-vs-chr17_$2-res_no-norm.txt

# Prr18-T vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:8530547:8636231 17 BP $2 $1_prr-t-vs-chr17_$2-res_no-norm.txt

# Pde10a vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:8977205:9077239 17 BP $2 $1_pde-vs-chr17_$2-res_no-norm.txt

# Qk vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:10462000:10562000 17 BP $2 $1_qk-vs-chr17_$2-res_no-norm.txt

# Zdhhc14 (extended to 100kb, start ~30kb upstream of CGI) vs chr17
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 17:5467100:5567180 17 BP $2 $1_zdh-vs-chr17_$2-res_no-norm.txt

# Xist (extended to 100kb) vs chrX
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic X:100615712:100718572 X BP $2 $1_xist-ext-vs-chrX_$2-res_no-norm.txt

# Kcnq1ot1 (extended 10kb at 5' end) vs chr7
#java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 7:150399016:150492452 7 BP $2 $1_kcn-vs-chr7_$2-res_no-norm.txt

# Rosa26 (extended to 100kb) vs chr6
java -jar ${directory}/juicer_tools.jar dump observed NONE $1.hic 6:113022310:113122310 6 BP $2 $1_rosa26-100kb-vs-chr6_$2-res_obs_no-norm.txt
java -jar ${directory}/juicer_tools.jar dump oe NONE $1.hic 6:113022310:113122310 6 BP $2 $1_rosa26-100kb-vs-chr6_$2-res_oe_no-norm.txt