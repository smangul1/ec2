#!/bin/bash


# file with sample commands for all the tools


FILE1=reads/input_1.fastq # first read file
FILE2=reads/input_2.fastq # second read file
KMER=18                   # kmer length
KMER2=3000000             # estimated number of all kmers
GLEN=3000                 # estimated genome length
RLEN=100                  # maximum read length


# create temporary directories

for tool in bfc bless coral fiona lighter musket pollux racer reckoner sga soapec; do 
    mkdir -p tmp/$tool; 
done



# run bfc
bash run.bfc.sh $FILE1 $FILE2 tmp/bfc $KMER $GLEN

# run bless
bash run.bless.sh $FILE1 $FILE2 tmp/bless/ $KMER

# run coral
bash run.coral.sh $FILE1 $FILE2 tmp/coral $KMER

# run fiona
bash run.fiona.sh $FILE1 $FILE2 tmp/fiona $GLEN

# run lighter
bash run.lighter.sh $FILE1 $FILE2 tmp/lighter/ $KMER $KMER2

# run musket
bash run.musket.sh $FILE1 $FILE2 tmp/musket/ $KMER $KMER2

# run pollux
bash run.pollux.sh $FILE1 $FILE2 tmp/pollux/ $KMER

# run racer
bash run.racer.sh $FILE1 $FILE2 tmp/racer/ $GLEN

# run reckoner
bash run.reckoner.sh $FILE1 $FILE2 tmp/reckoner/ $KMER $GLEN

# run sga
bash run.sga.sh $FILE1 $FILE2 tmp/sga/ $KMER

# run soapec
bash run.soapec.sh $FILE1 $FILE2 tmp/soapec/ $KMER $RLEN





