from multiprocessing import Pool
from tqdm import tqdm
from abnumber import Chain
from functools import partial
AHO_MAX_LEN = 152  # AHo has a fixed maximum length of sequences
from constants import PARENTAL_ID_TO_AA_SEQ


def get_cdrs(sequence, gap_character="-"):
    chain = Chain(sequence, scheme="aho", cdr_definition="chothia")

    cdr1_seq = chain.cdr1_seq
    cdr2_seq = chain.cdr2_seq
    cdr3_seq = chain.cdr3_seq

    cdr_indexes_dict = {
        1:(int(sequence.find(cdr1_seq)), int(sequence.find(cdr1_seq) + len(cdr1_seq))),
        2:(int(sequence.find(cdr2_seq)), int(sequence.find(cdr2_seq) + len(cdr2_seq))),
        3:(int(sequence.find(cdr3_seq)), int(sequence.find(cdr3_seq) + len(cdr3_seq))),
    }

    cdr_seqs_dict = {
        1:cdr1_seq,
        2:cdr2_seq,
        3:cdr3_seq,
    }
    
    for key in cdr_indexes_dict:
        start = cdr_indexes_dict[key][0]
        end = cdr_indexes_dict[key][1]
        assert cdr_seqs_dict[key] == sequence[cdr_indexes_dict[key][0]:cdr_indexes_dict[key][1]]

    return cdr_indexes_dict, cdr_seqs_dict


if __name__ == "__main__":
    parentals = ["A34467", "A38781" ]
    for parental in parentals:
        print(f"\n {parental}")
    
        seq = PARENTAL_ID_TO_AA_SEQ[parental]

        idxs, seqs = get_cdrs(sequence=seq, gap_character="-")
        for cdr_n in [1,2,3]:
            print("cdr ", cdr_n, seqs[cdr_n], idxs[cdr_n])

# A34467
# cdr  1 GRMFSAY (25, 32)
# cdr  2 TSGGN (51, 56)
# cdr  3 RDSVRYRDY (97, 106)

#  A38781
# cdr  1 GSVSDTN (25, 32)
# cdr  2 TTGGL (51, 56)
# cdr  3 RDIFRGNEY (97, 106)
