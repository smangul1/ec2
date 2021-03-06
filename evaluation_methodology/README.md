
# Statistical measures

We define the follwoing statistical measures

## Read level 

- TP - a successfully corrected read, i.e. all the erros were fixed
- FP - read was with no errors and we introduced errors in the read
- FN - read was not corrected, read was with errors before and remain with errors. Thsi may corresponds to the scenario when some errors were corrected, but no tall of them 
- TN - read was untouched, it was with no errors and remains woth no errors afetr error correction

When we define PPV, Sensitivity, Gain and Accuracy  

# How Statistical measures are computed across datasets

## Rep_Seq_sim

We simulate the TCR transcrips as described in here Mangul, Serghei, et al. "Profiling adaptive immune repertoires across multiple human tissues by RNA Sequencing." bioRxiv (2017): 089235.

In total we simulated 100 transcripts. We transcripts are available here

```
/u/home/i/imandric/project-zarlab/ErrorCorrection/DATASETS/S2_dataset/reference/[IGH|TRA]/trans.fa
```


We simulate reads from the IGH and TCRA transcripts (reffered as immune transcripts) using this command 

```
        ./simNGS/bin/simLibrary -r $l -i $ins -x $k imrep_revision2/$genename/trans.fa > imrep_revision2/$genename/ref.frag
        ./simNGS/bin/simNGS runfiles/$l.runfile imrep_revision2/$genename/ref.frag -p paired > imrep_revision2/$genename/sim_rl_${l}_cov_${k}.fq
```
Here -r is read length, -i insert size, -x coverage, genename is either TRA or IGH.


Simulator simulates reads with XXX error rate. More precisely, it uses Illumina profile files which originate from the sequencing machine. 
To obtain the error-free version of the reads (golden standard data), we map reads onto the immune transcripts using this command 

```
bowtie2-build -q ${transcriptFile} ${indexDir}/index
bowtie2 --quiet --no-hd --reorder -k 10 -q -p 10 -x ${indexDir}/index -U ${reads} -S reads.sam
```


To obtain the error-free version of the read we correct all mistahes with the reference using this script

```
/u/home/i/imandric/project-zarlab/ErrorCorrection/DATASETS/S2_dataset/ground_truth.py
```

After we obtain golden standard reads, we use the following scripts to compute the statistical measures 

```
https://github.com/smangul1/ec/tree/master/scripts
```

Namely, the scripts are: check_ED_bless.py,check_GT_bless.py, gain.py.

The script outputs 
Number of TP, FP, FN, TN and Gain and Accuracy


## WGS_sim

We simulate reads from chr1 of human genome (release 89), using this command 

```
cd /u/home/i/imandric/project-zarlab/ErrorCorrection/DATASETS/wgs_simulation/scripts_and_logs
```
```
../wgsim/wgsim -r 0.001 -R 0.0001 -e 0.005 -1 100 -2 100 -A 0 -N 1265528 ../reference/Homo_sapiens.GRCh38.dna.chromosome.1.fa ../datasets/t1/wgsim_rl_100_cov_1.1.fastq ../datasets/t1/wgsim_rl_100_cov_1.2.fastq > ../scripts_and_logs/t1/snps_rl_100_cov_1.txt
```
Note that the last command is in the file run.t1.sim_rl_100_1.sh in the same directory.
Explanation:
- a) -r 0.001 -> mutation rate 0.001
- b) -R 0.0001 -> fraction of indels
- c) -e 0.005 -> base error rate
- d) -1 100 -2 100 -> read length 100
- e) -N 1265528 - number of reads


We simulate 2 different datasets at different error rates. 
- T1 error rate: -r 0.001, -R 0.0001, -e 0.005
- T3  error rate: -r 0.03, -R 0.005, -e 0.02

It seems that I forgot to set indels fraction to 0. That's ok, I just avoid checking these reads.



The list of SNPs is available here

```
/u/home/i/imandric/project-zarlab/ErrorCorrection/DATASETS/wgs_simulation/scripts_and_logs/t1
```

and

```
/u/home/i/imandric/project-zarlab/ErrorCorrection/DATASETS/wgs_simulation/scripts_and_logs/t3
```



# RepSeq_real

We obtain TCR-Seq data from here

```
https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP045430
 ```


Given a file with UMIs we run this scripts to group the UMIs and correct the errors. Those reads are considered as golden standard reads. 

```
/u/home/n/ngcrawfo/project-zarlab/igor/VIRAL_CORRECTION/tools/data/datasets_no_barcodes_no_index/for_paper/ground_truth_TCRSEQ.py $READS_NO_BARCODE $BARCODE_FASTA_FILE $OUTPUT_FILE $CLUSTER_THRESHOLD $DISREGARD_FREQUENCY
```
Here CLUSTER_THRESHOLD is 5 by default (minimal nr of reads in a cluster) and DISREGARD_FREQUENCY is if the most popular nucleotide frequency less than this, drop the cluster.


The UMIs are here
```
/u/home/n/ngcrawfo/project-zarlab/igor/VIRAL_CORRECTION/tools/data/datasets_no_barcodes_no_index/for_paper/barcodes
```


The original reads without UMIs are here
```
/u/home/n/ngcrawfo/project-zarlab/igor/VIRAL_CORRECTION/tools/data/datasets_no_barcodes_no_index/for_paper/without_barcode
```






















