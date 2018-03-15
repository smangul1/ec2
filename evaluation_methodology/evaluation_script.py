
######################################################################
###Evaluation Code for Error Correction Benchmarking
#  project zar-lab ucla
#  3/14/18

#  Functions Contained: analyze_corrections,
######################################################################

#small testing code below the functions feel free to try to play around with it.


from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def analyze_corrections(true_sequence, raw_sequence, ec_sequence):
    """Keith Mitchell (keithgmitchell@g.ucla.edu).
        Function: Analyze true, raw(error_prone), and flawed/flawless error corrected reads
        Input 1: (string) true_sequence to be compared to raw and tool error corrected seq.
        Input 2: (string) raw_sequence (error prone sequence) that is to be compare to true and error corrected seq.
        Input 3: (string) ec_sequence (error corrected sequence) to be analyzed for flawed/flawless error correction
        Returns:
            Object 1: (dictionary) stats_dict which has the base level counts for the sequence:
                EX:        stats_dict = {'TN':0, 'TP':0, 'FN':0, 'FP':0, 'INDEL':0, 'TRIM': 0}
            Object 2: (string) seq_ID which dictates the read level classification (one of the following list):
                EX:        ["TN", "FP(NORMAL)", "FP(INDEL)", "FP(TRIMMING)", "TP", "FN"]
    """
    alignments = pairwise2.align.globalms(true_sequence, raw_sequence, 5, -4, -3, -0.1)

    true_aligned_raw = alignments[0][0]
    raw_aligned_true = alignments[0][1]

    alignments_2 = pairwise2.align.globalms(true_aligned_raw, ec_sequence, 5, -4, -3, -0.1)

    true_aligned_raw_ec = alignments_2[0][0]
    ec_aligned_true = alignments_2[0][1]

    alignments_3 = pairwise2.align.globalms(ec_aligned_true, raw_aligned_true, 5, -4, -3, -0.1)

    ec_aligned_final = alignments_2[0][0]
    raw_aligned_final = alignments_2[0][1]

    print "True:", true_aligned_raw_ec
    print "Raw: ", raw_aligned_final
    print "EC:  ", ec_aligned_final

    stats_dict = {'TN':0, 'TP':0, 'FN':0, 'FP':0, 'INDEL':0, 'TRIM': 0}

    if len(true_aligned_raw_ec) == len(ec_aligned_final) == len(raw_aligned_final):
        for true_bp, raw_bp, ec_bp in zip(true_aligned_raw_ec, raw_aligned_final, ec_aligned_final):
            if true_bp == raw_bp == ec_bp:
                stats_dict['TN'] += 1
            elif true_bp == raw_bp != ec_bp:
                stats_dict['FP'] += 1
            elif true_bp != raw_bp != ec_bp and ec_bp == true_bp:
                stats_dict['TP'] += 1
            elif true_bp != raw_bp == ec_bp:
                stats_dict['FN'] += 1
            if (true_bp == "-") or (raw_bp == "-") or (ec_bp == "-"):
                stats_dict['INDEL'] += 1


        ##this portion decides what to classify the sequence as a whole as.

        seq_classes = ["TN", "FP(NORMAL)", "FP(INDEL)", "FP(TRIMMING)", "TP", "FN"]
        seq_ID = ""
        if stats_dict['TN'] == len(true_aligned_raw_ec):
            seq_ID = seq_classes[0]
        elif stats_dict['FP'] != 0 and stats_dict['TP'] != 0 and stats_dict['INDEL'] == 0:
            seq_ID = seq_classes[1]
        elif stats_dict['FP'] != 0 and stats_dict['INDEL'] != 0:
            seq_ID = seq_classes[2]
        elif stats_dict['TP'] != 0 and stats_dict['FP'] == 0:
            seq_ID = seq_classes[4]

        #what to do about number 7, also need to incorporate the "FN"

        print "BASE LEVEL: ", stats_dict
        print "READ LEVEL: ", seq_ID

    else:
        print true_aligned_raw_ec, raw_aligned_final, ec_aligned_final

    return stats_dict, seq_ID

### Aggregate the results from each read into a data structure that can easily incorporate into excel.
### using the tool, dataset, and sequence ID.




true_sequence = 'ACCGGTCCA'
raw_sequence =  'ACCGGTCCA'
ec_sequence =   'CCGTGTCC'

analyze_corrections(true_sequence, raw_sequence, ec_sequence)




