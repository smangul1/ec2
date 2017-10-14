#!/bin/bash


# file with sample commands for all the tools


FILE1=reads/input_1.fastq # first read file
#FILE2=reads/input_2.fastq # second read file
KMER=18                   # kmer length
KMER2=3000000             # estimated number of all kmers
GLEN=3000                 # estimated genome length
#RLEN=100                  # maximum read length


# create temporary directories

for tool in bfc bless coral fiona lighter musket pollux racer reckoner sga soapec; do 
    mkdir -p $PWD/tmp/${tool}; 
done



# run bfc
bash run.bfc.se.sh $FILE1 $PWD/tmp/bfc $KMER

# run bless
bash run.bless.se.sh $FILE1 $PWD/tmp/bless/ $KMER

# run coral
bash run.coral.se.sh $FILE1 $PWD/tmp/coral $KMER

# run fiona
bash run.fiona.se.sh $FILE1 $PWD/tmp/fiona $GLEN

# run lighter
bash run.lighter.se.sh $FILE1 $PWD/tmp/lighter/ $KMER $GLEN

# run musket
bash run.musket.se.sh $FILE1 $PWD/tmp/musket/ $KMER

# run pollux
bash run.pollux.se.sh $FILE1 $PWD/tmp/pollux/ $KMER

# run racer
bash run.racer.se.sh $FILE1 $PWD/tmp/racer/ $GLEN

# run reckoner
bash run.reckoner.se.sh $FILE1 $PWD/tmp/reckoner/ $KMER

# run sga
bash run.sga.se.sh $FILE1 $PWD/tmp/sga/ $KMER

# run soapec
bash run.soapec.se.sh $FILE1 $PWD/tmp/soapec/ $KMER
