

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline


######################################################################
###Evaluation Script Methodology
######################################################################
# True and Raw Sequences are all the same to compare to for all of the EC reads from the seperate tools.
# This might be a list of files since there are multiple data sets to evaluate.
#     true_file = "file_path.fastq"[]
#     raw_file = "file_path.fastq"[]

# Given those are the same we just need to iterate through the ec files and compare to the corresponding raw/true file.
#     RUSSEL see what you can do here
#     need to be sure we know which tool each EC file is from
#     create temp. fasta file with the three sequences in it.

# Perform MSA using muscle to create the align object using the fasta file created in the previous step.

# Use MSA to analyze the TN, TP, FP, FN for each read.

# Aggregate the results from each read into a data structure that can easily incorporate into excel.

# (might be beneficial to make some of the above steps into definitions in the script.)

######################################################################


### True and Raw Sequences are all the same to compare to for all of the EC reads from the seperate tools.
### Given those are the same we just need to iterate through the ec files and compare to the corresponding raw/true file.

#for dataset in datasets:
    #for tool in tools:

tool_name = ""
true = SeqRecord(Seq("", generic_dna), id="Alpha")
raw = SeqRecord(Seq("", generic_dna), id="Beta")
ec = SeqRecord(Seq("", generic_dna), id="Gamma")

#compile these three sequences of interest into a temp. file called unaligned.fasta

#can add annotations to sequences with tool name/dataset name here potentially using seqrecord library object.




### Perform MSA using muscle to create the align object using the fasta file created in the previous step.

muscle_exe = r"C:\Program Files\Aligments\muscle3.8.31_i86win32.exe"
in_file = "..\unaligned.fasta"
out_file = "..\aligned.fasta"
muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
print(muscle_cline)
"C:\Program Files\Aligments\muscle3.8.31_i86win32.exe" -in "C:\My Documents\unaligned.fasta" -out "C:\My Documents\aligned.fasta"




### Use MSA to analyze the TN, TP, FP, FN for each read using the "aligned.fasta" file

TN, TP, FN, FP, INDEL = 0
#need to also factor in trimming here somehow

for true_bp, raw_bp, ec_bp in align[0], align[1], align[2]:
    if true_bp == raw_bp == ec_bp:
        TN +=1
    elif true_bp == raw_bp != ec_bp:
        FP +=1
    elif true_bp != raw_bp != ec_bp and ec_bp == true_bp:
        TP +=1
    elif true_bp != raw_bp == ec_bp:
        FN +=1
    elif true_bp or raw_bp or ec_bp == "-":
        INDEL += 1

seq_classes = ["TN", "FP(NORMAL)", "FP(INDEL)", "FP(TRIMMING)", "TP", "FN"]
seq_ID = ""

if TN == len(align):
    seq_ID = seq_classes[0]
elif FP != 0 and TP != 0 and INDEL == 0:
    seq_ID = seq_classes[1]
elif FP != 0 and TP != 0 and INDEL != 0:
    seq_ID = seq_classes[2]
elif TP != 0 and FP == 0:
    seq_ID = seq_classes[4]
#what to do about number 7, also need to incorporate the "FN"






### Aggregate the results from each read into a data structure that can easily incorporate into excel.



