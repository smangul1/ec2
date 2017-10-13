#!/bin/bash


# file with sample commands for all the tools


FILE1=reads/input.fastq   # read file
FILE1=/u/home/s/serghei/project/EC_survey/RepSeq_sim/IGH_sim_rl_100_cov_1_R1.fastq   # read file
KMER=18                   # kmer length
KMER2=3000000             # estimated number of all kmers
GLEN=3000                 # estimated genome length
RLEN=100                  # maximum read length
tmp=/tmp/blhill

# create temporary directories

for tool in bfc bless coral fiona lighter musket pollux racer reckoner sga soapec; do 
    mkdir -p $tmp/$tool; 
done



# run bfc
bash run.bfc.se.sh $FILE1 $tmp/bfc $KMER $GLEN

# run bless
bash run.bless.se.sh $FILE1 $tmp/bless/ $KMER

# run coral
bash run.coral.se.sh $FILE1 $tmp/coral $KMER

# run fiona
bash run.fiona.se.sh $FILE1 $tmp/fiona $GLEN

# run lighter
bash run.lighter.se.sh $FILE1 $tmp/lighter/ $KMER $KMER2

# run musket
bash run.musket.se.sh $FILE1 $tmp/musket/ $KMER $KMER2

# run pollux
bash run.pollux.se.sh $FILE1 $tmp/pollux/ $KMER

# run racer
bash run.racer.se.sh $FILE1 $tmp/racer/ $GLEN

# run reckoner
bash run.reckoner.se.sh $FILE1 $tmp/reckoner/ $KMER $GLEN

# run sga
bash run.sga.se.sh $FILE1 $tmp/sga/ $KMER

# run soapec
bash run.soapec.se.sh $FILE1 $tmp/soapec/ $KMER $RLEN





