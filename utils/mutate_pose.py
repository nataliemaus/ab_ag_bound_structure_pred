import os 
import sys 
sys.path.append("../")
from constants import PARENTAL_ID_TO_AG_SEQ 
from utils.load_pose import load_pose 
import pyrosetta 
pyrosetta.init()


def mutate_residues(
    path_to_pdb,
    parental,
    mutant_positions=[1,2], 
    mutant_aas=["A","A"], 
    save_mutated_pdb=False,
):
    pose = load_pose(path_to_pdb=path_to_pdb,)
    ag_seq = PARENTAL_ID_TO_AG_SEQ[parental]

    # First, we check if the antigen comes first, and if so modify the mutation positions accordingly
    input_seq = pose.sequence()
    if input_seq.startswith(ag_seq):
        temp = []
        for og_mutant_pos in mutant_positions:
            temp.append(og_mutant_pos + len(ag_seq) )
        mutant_positions = temp  

    # make each mutaiton passed in 
    for ix, pos in enumerate(mutant_positions):
        pyrosetta.toolbox.mutants.mutate_residue(pose, pos, mutant_aas[ix])

    if save_mutated_pdb: 
        # use a filename indiciative of what mutations were made and where 
        mutant_aa_str = ""
        for mut_aa in mutant_aas: mutant_aa_str += mut_aa 
        mutant_pos_str = ""
        for mut_pos in mutant_positions: 
            mutant_pos_str += "_"
            mutant_pos_str += str(mut_pos) 
        save_mutated_pdb_path = path_to_pdb.replace(".pdb", f"_muts_{mutant_aa_str}_at{mutant_pos_str}.pdb") 
        # save pose with new mutations to pdb 
        pose.dump_pdb(save_mutated_pdb_path)
        # return mutated pose and path where file was saved 
        return pose, save_mutated_pdb_path

    return pose 

