
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

In total we simulated ? transcripts. We transcripts are available here


We simulate reads from the IGH and TCRA transcripts (reffered as immune transcripts) using this command 

```
???
```

Simulator simulates reads with ? error rate. To obtain the error-free version of the reads (golden standard data), we map reads onto the immune transcripts using this command 

```
???
```


To obtain the error-free version of the read we correct all mistahes with the reference using this script

```
???
```

After we obtain golden standard reads, we use the following scripts to compute the statistical measures 

```
????
```

The script outputs 
Number of TP, FP, FN, TN and Gain and Accuracy


## WGS_sim

We simulate reads from chr1 of human genome (release ?), using this command 

```
???
```

We simulate 2 different datasets at different error rates. 
- T1 error rate ?
- T3  error rate

The simulateror besides errros introduse SNPs. To obtain the error-free version of the reads (golden standard data), we map reads onto the refference genome using this command

```
???
```




















