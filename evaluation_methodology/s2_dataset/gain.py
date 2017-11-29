

##  Script for computing gain metric on reads
##  This is how to run it:
##      cd /u/home/i/imandric/project-zarlab/ErrorCorrection/accuracy_scripts
##      LD_LIBRARY_PATH=/u/home/i/imandric/project-zarlab/tools/Python-2.7.5 python gain2.py musket_IGH_sim_rl_100_cov_1_R1_19.corrected.fastq ../DATASETS/S2_dataset/reads/IGH/sim_rl_100_cov_1/input_1.fastq ../DATASETS/S2_dataset/reads/IGH/sim_rl_100_cov_1/input_2.fastq ../DATASETS/S2_dataset/alignments/sam/IGH/sim_rl_100_cov_1/input_1.sam ../DATASETS/S2_dataset/alignments/sam/IGH/sim_rl_100_cov_1/input_2.sam ../DATASETS/S2_dataset/reference/IGH/trans.fa
##  Example output:
##      musket_IGH_sim_rl_100_cov_1_R1_19.corrected.fastq -0.450568678915 0.672790901137 769 1284 374 1143 1658
##  Here, the gain is -0.45


import sys, os
import random, string
import re
from Bio import SeqIO, AlignIO
from Bio import pairwise2
#from Bio.Emboss import Applications as apps

file1 = sys.argv[1] # this is the corrected reads file1
file2 = sys.argv[2] # this is the original reads file1
file3 = sys.argv[3] # this is the original reads file2
samfile1 = sys.argv[4] # this is the alignment of corrected reads to the reference transcripts
samfile2 = sys.argv[5]
tranfile = sys.argv[6]



def reverse_complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases




# original reads
origreads = {}


for record in SeqIO.parse(file2, "fastq"):
    origreads[record.id] = str(record.seq)
   

for record in SeqIO.parse(file3, "fastq"):
    origreads[record.id] = str(record.seq)
#------------------------------------------




# reading reference transcripts
transcripts = {}
for record in SeqIO.parse(tranfile, "fasta"):
    transcripts[record.id] = str(record.seq)

#-----------------------------------------------

with open(samfile1) as f:
    s1 = f.readlines()
s1 = map(lambda x: x.strip().split(), s1)

with open(samfile2) as f:
    s2 = f.readlines()
s2 = map(lambda x: x.strip().split(), s2)


a = map(lambda z: [z[0], z[2], int(z[3]), z[5], z[9], z[1]], s1 + s2)


em = set()
genomic_regions = dict() # this is the dictionary containing genomic regions where each of the reads map


for line in a:

    if line[5] == "4":
        continue

    start = line[2]
    segments = map(int, re.split("[M|N|I|D|S]", line[3])[:-1])
    letters = [x for x in line[3] if x in ["M", "N", "I", "D", "S"]]

    intervals = []
    ref_coord = int(start)
    read_coord = 0
    for l, s in zip(letters, segments):
        if l == "S":
            intervals.append(("S", (read_coord, read_coord + s), (ref_coord, ref_coord)))
            read_coord += s
        elif l == "I":
            intervals.append(("I", (read_coord, read_coord + s), (ref_coord, ref_coord)))
            read_coord += s
        elif l == "D":
            intervals.append(("D", (read_coord, read_coord), (ref_coord, ref_coord + s)))
            ref_coord += s
        elif l == "M":
            intervals.append(("M", (read_coord, read_coord + s), (ref_coord, ref_coord + s)))
            read_coord += s
            ref_coord += s
        #elif l == "N": # this case is for spliced alignments
        #    intervals.append(("N", (read_coord, read_coord), (ref_coord, ref_coord + s)))
        #    ref_coord += s


    for letter, read_coords, ref_coords in intervals:

        if letter == "M":
            for u, v in zip(range(*read_coords), range(*ref_coords)):
                litera = transcripts[line[1]][v - 1]
                if line[4][u] != litera:
                    if line[5] in ["0", "256"]:
                        em.add("%s:%s:M" % (line[0], u))
                    else:
                        em.add("%s:%s:M" % (line[0], len(line[4]) - 1 - u))
        elif letter == "I":
            if line[5] in ["0", "256"]:
                em.add("%s:%s-%s:I" % (line[0], read_coords[0], read_coords[1]))
            else:
                em.add("%s:%s-%s:I" % (line[0], len(line[4]) - 1 - read_coords[1], len(line[4]) - 1 - read_coords[0]))
        elif letter == "D":
            if line[5] in ["0", "256"]:
                em.add("%s:%s-%s:D" % (line[0], read_coords[0], read_coords[1]))
            else:
                em.add("%s:%s-%s:D" % (line[0], len(line[4]) - 1 - read_coords[1], len(line[4]) - 1 - read_coords[0]))
        elif letter == "S":
            if line[5] in ["0", "256"]:
                em.add("%s:%s-%s:S" % (line[0], read_coords[0], read_coords[1]))
            else:
                em.add("%s:%s-%s:S" % (line[0], len(line[4]) - 1 - read_coords[1], len(line[4]) - 1 - read_coords[0]))

    genomic_region_start = intervals[0][2][0] - 1
    genomic_region_end = intervals[-1][2][-1] - 1
    genomic_region_seq = transcripts[line[1]][genomic_region_start:genomic_region_end]
    if line[5] in ["16", "272"]:
        genomic_region_seq = reverse_complement(genomic_region_seq)

    genomic_regions[line[0]] = genomic_region_seq


corrdict = {}
for record in SeqIO.parse(file1, "fastq"):
    corrdict[record.id] = str(record.seq)


common = set(genomic_regions.keys()) & set(corrdict.keys())

# format for errors
# mismatches: seq:position
# indels: seq:start-end

ec = set()

for key in common:
    alignment = pairwise2.align.globalms(corrdict[key], genomic_regions[key], 5, -1, -10, -2)[0]
    #print key, origreads[key], corrdict[key], genomic_regions[key], alignment
    allen = len(alignment[0])
    read_position = -1
    insertion = False
    deletion = False
    instart = None
    delstart = None
    for i in range(allen):
        if alignment[0][i] != "-":
            read_position += 1 # this position belongs to the corrected read
        if alignment[0][i] == alignment[1][i]:
            if insertion == True:
                ec.add("%s:%s-%s:I" % (key, instart, read_position))
                insertion = False
                instart = None
            if deletion == True:
                ec.add("%s:%s-%s:D" % (key, delstart, read_position))
                deletion = False
                delstart = None
        else: # here is an error
            if alignment[0][i] == "-": # deletion
                deletion = True
                if delstart is None:
                    delstart = read_position
            elif alignment[1][i] == "-": #insertion
                insertion = True
                if instart is None:
                    instart = read_position
            else: #mismatch
                ec.add("%s:%s:M" % (key, read_position))

#for x in em:
#    print "EM", x


#for x in ec:
#    print "EC", x


tp = len(em - ec)
fp = len(ec - em)
fn = len(em & ec)

gain = (tp - fp) * 1.0 / (tp + fn)

sens = tp * 1.0 / (tp + fn)


print file1, gain, sens, tp, fp, fn, len(em), len(ec)

