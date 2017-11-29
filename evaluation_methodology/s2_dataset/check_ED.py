


## Script for computing edit distance histrogram 
## Input: ground truth file (.gt), corrected fastq file
## Output: edit distance histogram
## Example run: 
##     # cd /u/home/i/imandric/project-zarlab/ErrorCorrection/accuracy_scripts
##     # python check_ED.py ../DATASETS/S2_dataset/ground_truths/IGH/sim_rl_100_cov_1.gt  musket_IGH_sim_rl_100_cov_1_R1_19.corrected.fastq
## Output:
##     # 0:659 1:180 2:106 3:90 4:64 5:32 6:20 7:24 8:12 9:5 10:7 11:3 12:3 14:1 16:1 17:1



from Bio import SeqIO
import sys
import jellyfish


gt1 = sys.argv[1]
corrected_reads = sys.argv[2]


with open(gt1) as f:
    gt1 = f.readlines()

gt1 = map(lambda x: x.strip().split(), gt1)

gt1dict = {}

for x, y, z in gt1:
    gt1dict[x] = (y, z)

histogram = {}

rw = {}
handle = open(corrected_reads, "rU")
for record in SeqIO.parse(handle, "fastq"):
    if True:
        if record.id in gt1dict:
            corread = record.seq.tostring()
            origread, trueread = gt1dict[record.id]
            distance = jellyfish.levenshtein_distance(unicode(corread), unicode(trueread))
            if distance not in histogram:
                histogram[distance] = 0
            histogram[distance] += 1
handle.close()

hist = []
for key in sorted(histogram.keys()):
    hist.append("%s:%s" % (key, histogram[key]))

print " ".join(hist)

