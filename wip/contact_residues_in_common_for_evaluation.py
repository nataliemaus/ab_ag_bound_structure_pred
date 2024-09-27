from Bio import PDB
import glob 
import sys 
sys.path.append("../")
from constants import PARENTAL_ID_TO_AB_CHAIN_ID, PARENTAL_ID_TO_ANTIGEN_CHAIN_ID
import pandas as pd

def get_all_residue_pair_distances(
    input_pdb_filename, 
    ab_chain_id,
    ag_chain_id,
):
    # Parse the input PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb_filename)

    interface_residues = set() 

    # rn1, rn2, distance 
    df = {
        f"res_num_{ab_chain_id}":[],
        f"res_num_{ag_chain_id}":[],
        "distance":[],
    }

    analysis_done = False 
    for model in structure:
        for chain1 in model:
            for chain2 in model:
                if (chain1 != chain2) and (not analysis_done):
                    for residue1 in chain1:
                        res_num1 = get_res_num(residue_object=residue1)
                        chain_id1 = residue1.parent.id 
                        for residue2 in chain2:
                            res_num2 = get_res_num(residue_object=residue2)
                            chain_id2 = residue2.parent.id 
                            # Calculate minimum distance between any atoms of residue1 and residue2
                            min_distance = min(atom1 - atom2 for atom1 in residue1 for atom2 in residue2)
                            df[f"res_num_{chain_id1}"].append(res_num1)
                            df[f"res_num_{chain_id2}"].append(res_num2)
                            df["distance"].append(min_distance)
                    analysis_done = True 

    df = pd.DataFrame.from_dict(df)    # (13189, 5)
    df = df.sort_values(by=['distance']) # (13189, 5)

    csv_filename = input_pdb_filename.split("/")[-1].replace(".pdb", f"_ag{ag_chain_id}_ab{ab_chain_id}.csv")
    save_csv_path = f"contact_residues_csvs/{csv_filename}"
    df.to_csv(save_csv_path, index=False)

    return df



def get_res_num(residue_object):
    # get residue number 
    residue_info = residue_object.__repr__()
    res_num = None 
    for str1 in residue_info.split(" "):
        if "resseq=" in str1: 
            res_num = str1.replace("resseq=", "")
            res_num = int(res_num)
            break
    if res_num is None: 
        assert 0, "error grabbing residue number"
    return res_num 


# TODO: 
# 1. take in list of pdb file paths 
# 2. compute contact residues on ab AND ag 
# 3. compute contact residues common amont all files, return 
# 4. use those contact residues in common to pass into HDock! 

def save_all_residue_pair_distances(
    known_pose_path,
    pred_pose_path,
    ab_chain_id="A",
    ag_chain_id="B",
    known_dist_threshold=5.0,
    pred_dist_threshold=8.0,
):
    df_known = get_all_residue_pair_distances(
        input_pdb_filename=known_pose_path,
        ab_chain_id=ab_chain_id,
        ag_chain_id=ag_chain_id,
    )
    df_known = df_known[df_known["distance"] <= known_dist_threshold]
    df_pred = get_all_residue_pair_distances(
        input_pdb_filename=pred_pose_path,
        ab_chain_id=ab_chain_id,
        ag_chain_id=ag_chain_id,
    )
    df_pred = df_pred[df_pred["distance"] <= pred_dist_threshold]

    known_ab_contact_residues = df_known[f"res_num_{ab_chain_id}"].values.tolist()
    known_ag_contact_residues = df_known[f"res_num_{ag_chain_id}"].values.tolist()

    pred_ab_contact_residues = df_pred[f"res_num_{ab_chain_id}"].values.tolist()
    pred_ag_contact_residues = df_pred[f"res_num_{ag_chain_id}"].values.tolist()
    
    n_good = 0
    n_total = 0
    for res_num in known_ab_contact_residues:
        if res_num in pred_ab_contact_residues:
            n_good += 1 
        n_total += 1
    
    for res_num in known_ag_contact_residues:
        if res_num in pred_ag_contact_residues:
            n_good += 1 
        n_total += 1

    perc_good = n_good/n_total
    print(f"    {n_good}/{n_total} = {perc_good} Contact Residues in Known Pose are also in Pred")


if __name__ == "__main__":
    parental_to_best_round2_hdock_model_num = {
        "A34467":34, 
        "A38781":37, 
    }
    
    for parental in parental_to_best_round2_hdock_model_num.keys():
        print(f"PARENTAL {parental}::: ")
        known_pose_path = f"pdbs/{parental}_all_mutations_known_pose/maus_{parental}_known_pose.pdb" 
        model_num = parental_to_best_round2_hdock_model_num[parental]
        pred_pose_path = f"pdbs/hdock_round2_{parental}_binding_sites_10_0_chainids_modified/model_{model_num}_chain_ids_modified.pdb" 
        save_all_residue_pair_distances(
            known_pose_path=known_pose_path,
            pred_pose_path=pred_pose_path,
            ab_chain_id=PARENTAL_ID_TO_AB_CHAIN_ID[parental], 
            ag_chain_id=PARENTAL_ID_TO_ANTIGEN_CHAIN_ID[parental], 
            known_dist_threshold=5.0,
            pred_dist_threshold=8.0,
        )


# PARENTAL A34467::: 
#     93/156 = 0.5961538461538461 Contact Residues in Known Pose are also in Pred
# PARENTAL A38781::: 
#     128/134 = 0.9552238805970149 Contact Residues in Known Pose are also in Pred