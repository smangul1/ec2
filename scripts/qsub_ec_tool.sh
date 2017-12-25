#!/bin/bash

tool="sga"
script_dir=$HOME/code/github/ec2/wrappers/paired_end
input_file_list=$HOME/code/github/ec2/datasets/table_data_PE.csv
output_dir=$SCRATCH/ec/${tool}

#for k in 15 19 24 27 31
# for each k-mer 
for k in 19 
do
    # every other line in file is a separate dataset
    for i in {2..206..2}
    do
        echo "k: $k i: $i"
        ls $output_dir/k${k}_$i/*.corrected.fastq.gz
        if [ "$?" -eq "0" ]; then
            echo "corrected file exists... skipping"
            continue
        fi
        if_1=`sed "${i}q;d" $input_file_list | cut -f1 -d ','`
        if_2=`sed "$((i+1))q;d" $input_file_list | cut -f1 -d ','`
        echo "$if_1 $if_2"
        qsub -V -N "${tool}_$k_$i" -l h_data=16G,highp,time=03:59:00 $script_dir/run.${tool}.sh $if_1 $if_2 "$output_dir/k${k}_$i" $k
    done
done
