import sys 
sys.path.append("../")
from utils.mutate_utils import get_mutations_convert_seq1_to_seq2

def get_mutations(
    parental_seq,
    seqs_list,
):
    mutatant_positions_per_seq, mutatant_aas_per_seq = [], []
    for mutated_seq in seqs_list:
        positions_list, aas_list, _ = get_mutations_convert_seq1_to_seq2(
            seq1=parental_seq,
            seq2=mutated_seq,
        )
        mutatant_positions_per_seq.append(positions_list)
        mutatant_aas_per_seq.append(aas_list)

    return mutatant_positions_per_seq, mutatant_aas_per_seq 