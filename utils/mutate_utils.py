from Bio import pairwise2

def get_mutations_convert_seq1_to_seq2(
    seq1,
    seq2,
):
    # i.e. convert parental (seq1) to mutated seq (seq2)
    positions_list = []
    aas_list = []
    alignment = pairwise2.align.globalxx(seq1, seq2)[0]

    aligned_seq1 = alignment.seqA
    aligned_seq2 = alignment.seqB 
    n_mutations1 = aligned_seq1.count("-") # parental
    n_mutations2 = aligned_seq2.count("-")
    if n_mutations1 != n_mutations2:
        # must be insertions or deletions, not just mutations, ignore 
        return [], [], 1_000 
    for ix, char in enumerate(aligned_seq1):
        if char == "-":
            positions_list.append(ix + 1 - n_mutations1)
            aas_list.append(aligned_seq2[ix])
    # positions_list: list of integers giving postiions where we mutate seq1 to get seq2
    # aas_list: list of chars giving new aa's to mutate in to end up with seq2 
    n_mutations = len(aas_list)
    return positions_list, aas_list, n_mutations 