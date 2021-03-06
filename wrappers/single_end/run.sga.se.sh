#!/bin/bash

AUTHOR="imandric1"



################################################################
##########          The main template script          ##########
################################################################

toolName="sga"
toolPath="/u/home/b/blhill/local/bin/sga"

# STEPS OF THE SCRIPT
# 1) prepare input if necessary
# 2) run the tool
# 3) transform output if necessary
# 4) compress output


# THE COMMAND LINE INTERFACE OF THE WRAPPER SCRIPT
# $tool $input1 $input2 $outdir $kmers $others
# |      mandatory part       | | extra part |
# <---------------------------> <------------>




if [ $# -lt 3 ]
then
echo "********************************************************************"
echo "Script was written for project : Best practices for conducting benchmarking in the most comprehensive and reproducible way"
echo "This script was written by Igor Mandric"
echo "********************************************************************"
echo ""
echo "1 <input>  - .fastq"
echo "2 <outdir> - dir to save the output"
echo "3 <kmer>   - kmer length"
echo "--------------------------------------"
exit 1
fi



# mandatory part
input=$1
outdir=$2

# extra part (tool specific)
kmer=$3



# STEP 0 - create output directory if it does not exist

mkdir -p $outdir

# -----------------------------------------------------


# STEP 1 - prepare input if necessary (ATTENTION: TOOL SPECIFIC PART!)

pwd=$PWD
cd $outdir
outdir_abs=$PWD
cd $pwd

logfile=$outdir/report_$(basename ${input%.*})_${toolName}_${kmer}.log

echo "START" > $logfile

# -----------------------------------




# STEP 2 - run the tool (ATTENTION: TOOL SPECIFIC PART!)

now="$(date)"
printf "%s --- RUNNING %s\n" "$now" $toolName >> $logfile

# run the command
res1=$(date +%s.%N)

pwd=$PWD
cd $outdir_abs
$toolPath preprocess $input -o one_output_file.preprocessed.fastq >> $logfile 2>&1
$toolPath index -a ropebwt one_output_file.preprocessed.fastq >> $logfile 2>&1
$toolPath correct -k $kmer -o one_output_file.out.fastq one_output_file.preprocessed.fastq >> $logfile 2>&1


res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
now="$(date)"
printf "%s --- TOTAL RUNTIME: %d:%02d:%02d:%02.4f\n" "$now" $dd $dh $dm $ds >> $logfile

now="$(date)"
printf "%s --- FINISHED RUNNING %s %s\n" "$now" $toolName >> $logfile

# ---------------------




# STEP 3 - transform output if necessary (ATTENTION: TOOL SPECIFIC PART!)

# if you need to transform fasta to fastq - here is the command:
#     awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' file.fastq > file.fasta

now="$(date)"
printf "%s --- TRANSFORMING OUTPUT\n" "$now" >> $logfile


cat one_output_file.out.fastq | gzip > ${toolName}_$(basename ${input%.*}).corrected.fastq.gz
rm one_output_file.preprocessed.fastq
rm one_output_file.out.fastq
rm one_output_file.*

now="$(date)"
printf "%s --- TRANSFORMING OUTPUT DONE\n" "$now" >> $logfile

# --------------------------------------
cd $pwd


printf "DONE" >> $logfile



