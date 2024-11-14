import pandas as pd 
import glob 
import numpy as np 
import os 
import copy 
from utils.get_mutations import get_mutations
from utils.refine_pose import refine_pose 
from utils.mutate_pose import mutate_residues 
from utils.score_pose import get_score_function
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot
from constants import PARENTAL_ID_TO_AB_SEQ, HER2_CDR3_FIRST_INDEX, HER2_CDR3_LAST_INDEX

# Deubg --> do fewer sequences and poses to more quickly check the the whole pipeline works 
DEBUG_MODE = True 


def main(
    paths_to_poses_list,
    parental,
    affinity_data_path,
    parental_seq,
    save_correlations_csv_path,
    save_plots_dir,
    make_scatter_plots=False,
    skip_refinement=False,
):
    df = pd.read_csv(affinity_data_path)
    # df = df[df["n_mutations"] > 0] # remove first row with 0 mutations (un-mutated parental)

    if DEBUG_MODE:
        df = df[0:4]

    # kd_per_seq = df["mean_neg_log_kd"].values.tolist() # --> affinity_per_seq
    # seqs_list = df["seq"].values.tolist()

    new_cdr3s = df["seq"].values.tolist()
    seqs_list = []
    for new_cdr3 in new_cdr3s:
        new_seq = copy.deepcopy(parental_seq)
        new_seq = [char for char in new_seq]
        new_cdr3 = [char for char in new_cdr3]
        new_seq[HER2_CDR3_FIRST_INDEX:HER2_CDR3_LAST_INDEX] = new_cdr3 
        new_seq = "".join(new_seq)
        seqs_list.append(new_seq)
    

    affinity_per_seq = df["label"].values.astype(np.float32) 
    affinity_per_seq[df["class"].values == "mid"] = 0.5 # high --> 1, low --> 0, mid --> 0.5 
    affinity_per_seq = affinity_per_seq.tolist()

    mutatant_positions_per_seq, mutatant_aas_per_seq = get_mutations(
        parental_seq=parental_seq,
        seqs_list=seqs_list,
    )

    energy_sfxn = get_score_function()
    # for each pose in list, compute kd-energy spearman correlation coeff 
    correlation_coeffs = []
    for pose_path in paths_to_poses_list:
        binding_energies_per_seq = []
        for ix, mutatnt_positions_list in enumerate(mutatant_positions_per_seq):
            mutated_pose = mutate_residues(
                path_to_pdb=pose_path,
                parental=parental,
                mutant_positions=mutatnt_positions_list,
                mutant_aas=mutatant_aas_per_seq[ix],
                save_mutated_pdb=False,
            )
            if DEBUG_MODE:
                print("mutated pose:", mutated_pose)
            if skip_refinement:
                refined_pose = mutated_pose
            else:
                refined_pose = refine_pose(
                    parental=parental,
                    pose=mutated_pose,
                ) 
            if DEBUG_MODE:
                print("refined pose:", refined_pose)
            binding_energy = energy_sfxn(refined_pose) 
            if DEBUG_MODE:
                print("binding_energy:", binding_energy)
            binding_energies_per_seq.append(binding_energy)
        binding_energies_per_seq = np.array(binding_energies_per_seq)
        spearman_r = get_spearman_r(affinity_per_seq,binding_energies_per_seq)
        if DEBUG_MODE:
            print("final binding_energies_per_seq:", binding_energies_per_seq)
            print("spearman_r:", spearman_r) 
        if make_scatter_plots:
            plot_filename = pose_path.split("/")[-1].replace(".pdb", ".png")
            rounded_spearman_r = round(spearman_r, 3)
            if not os.path.exists(save_plots_dir):
                os.mkdir(save_plots_dir) 
            scatter_plot( 
                x=binding_energies_per_seq, 
                y=affinity_per_seq, 
                x_label="Rosetta Binding Energy", 
                y_label="Affinity", 
                title=f"KD-Energy Correlation Spearman Coeff = {rounded_spearman_r}", 
                path_to_save_plot=f"{save_plots_dir}{plot_filename}", 
            )
        correlation_coeffs.append(spearman_r)
    

    save_df = {}
    save_df["kd_energy_spearman_r"] = np.array(correlation_coeffs)
    save_df["pose_path"] = np.array(paths_to_poses_list)
    save_df = pd.DataFrame.from_dict(save_df)
    save_df = save_df.sort_values(by=["kd_energy_spearman_r"], ascending=False)
    save_df.to_csv(save_correlations_csv_path, index=None)


if __name__ == "__main__":
    parental = "HER2"
    skip_refinement = False 
    if skip_refinement:
        results_dir = "her2/results_no_refinement"
    else:
        results_dir = "her2/results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir) 
    affinity_data_path = f"her2/affinity_data_cdr3.csv"
    parental_seq = PARENTAL_ID_TO_AB_SEQ[parental]
    save_correlations_csv_path = f"{results_dir}/{parental}_hdock_round1_poses_correlation_results.csv"
    paths_to_round1_hdock_predicted_poses = glob.glob(f"her2/hdock_her2/model_*.pdb")
    # remove .clean.pdb versions 
    temp = []
    for hdock_pdb_file in paths_to_round1_hdock_predicted_poses:
        if not ("clean" in hdock_pdb_file):
            temp.append(hdock_pdb_file)
    paths_to_round1_hdock_predicted_poses = temp 
    if DEBUG_MODE:
        paths_to_round1_hdock_predicted_poses = paths_to_round1_hdock_predicted_poses[0:3]
    main(
        paths_to_poses_list=paths_to_round1_hdock_predicted_poses,
        parental=parental,
        affinity_data_path=affinity_data_path,
        parental_seq=parental_seq,
        save_correlations_csv_path=save_correlations_csv_path,
        make_scatter_plots=True,
        save_plots_dir=f"{results_dir}/plots/",
        skip_refinement=skip_refinement,
    )

    # Ran in debug mode with refinement to confirm this works now that I'm not passing in lig.pdb file :) yay
        # tmux new -s pillbox
        # python3 run_pipeline.py 
        # tmux new -s her2
        # python3 run_pipeline_her2.py   (debug mode running)