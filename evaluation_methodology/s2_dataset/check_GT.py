

## Script for computing TP, FP, etc 
## Input: ground truth file (.gt), corrected fastq file
## Output: read level metrics
## Example run: 
##     # cd /u/home/i/imandric/project-zarlab/ErrorCorrection/accuracy_scripts
##     # python check_GT.py ../DATASETS/S2_dataset/ground_truths/IGH/sim_rl_100_cov_1.gt  musket_IGH_sim_rl_100_cov_1_R1_19.corrected.fastq
## Output:
##     # corrected: 240 0.198675496689
##     # remained: 419 0.346854304636
##     # undercorrected: 275 0.227649006623
##     # overcorrected: 274 0.226821192053


from Bio import SeqIO
import sys


gt1 = sys.argv[1]
corrected_reads = sys.argv[2]


with open(gt1) as f:
    gt1 = f.readlines()


gt1 = map(lambda x: x.strip().split(), gt1)

gt1dict = {}

for x, y, z in gt1:
    gt1dict[x] = (y, z)


corrected = 0
remained = 0
undercorrected = 0
overcorrected = 0

rw = {}
handle = open(corrected_reads, "rU")
for record in SeqIO.parse(handle, "fastq"):
    if True:
        if record.id in gt1dict:
            corread = record.seq.tostring()
            origread, trueread = gt1dict[record.id]
            if corread == trueread:
                if origread == trueread:
                    remained += 1
                else:
                    corrected += 1
            else:
                if trueread == origread:
                    overcorrected += 1
                else:
                    undercorrected += 1

handle.close()

sumall = float(corrected + remained + undercorrected + overcorrected)

print "corrected:", corrected, corrected / sumall
print "remained:", remained, remained / sumall
print "undercorrected:", undercorrected, undercorrected / sumall
print "overcorrected:", overcorrected, overcorrected / sumall

