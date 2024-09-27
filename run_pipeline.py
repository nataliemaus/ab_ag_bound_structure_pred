import pandas as pd 
import glob 
import numpy as np 
import os 
from utils.mutate_utils import get_mutations_convert_seq1_to_seq2
from utils.refine_pose import refine_pose 
from utils.mutate_pose import mutate_residues 
from utils.score_pose import get_score_function
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot
from constants import PARENTAL_ID_TO_AA_SEQ

# Deubg --> do fewer sequences and poses to more quickly check the the whole pipeline works 
DEBUG_MODE = False 

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


def main(
    paths_to_poses_list,
    cfps_data_path,
    parental_seq,
    save_correlations_csv_path,
    save_plots_dir,
    make_scatter_plots=False,
):
    df = pd.read_csv(cfps_data_path)
    df = df[df["n_mutations"] > 0] # remove first row with 0 mutations (un-mutated parental)

    if DEBUG_MODE:
        df = df[0:10]

    kd_per_seq = df["mean_neg_log_kd"].values.tolist()
    seqs_list = df["seq"].values.tolist()

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
                mutant_positions=mutatnt_positions_list,
                mutant_aas=mutatant_aas_per_seq[ix],
                save_mutated_pdb=False,
            )
            refined_pose = refine_pose(
                pose=mutated_pose,
            ) 
            binding_energy = energy_sfxn(refined_pose) 
            binding_energies_per_seq.append(binding_energy)
        binding_energies_per_seq = np.array(binding_energies_per_seq)
        spearman_r = get_spearman_r(kd_per_seq,binding_energies_per_seq)
        if make_scatter_plots:
            plot_filename = pose_path.split("/")[-1].replace(".pdb", ".png")
            rounded_spearman_r = round(spearman_r, 3)
            if not os.path.exists(save_plots_dir):
                os.mkdir(save_plots_dir) 
            scatter_plot( 
                x=binding_energies_per_seq, 
                y=kd_per_seq, 
                x_label="Rosetta Binding Energy", 
                y_label="Negative Log KD (CFPS)", 
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
    results_dir = "pillbox/results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir) 
    relevant_parentals = ["A38781", "A34467"]
    for parental in relevant_parentals:
        cfps_data_path = f"pillbox/pillbox_{parental}_human_cfps.csv"
        parental_seq = PARENTAL_ID_TO_AA_SEQ[parental]
        parental_results_dir = f"{results_dir}/{parental}"
        if not os.path.exists(parental_results_dir):
            os.mkdir(parental_results_dir) 
        save_correlations_csv_path = f"{parental_results_dir}/{parental}_hdock_round1_poses_correlation_results.csv"
        paths_to_round1_hdock_predicted_poses = glob.glob(f"pillbox/HDOCK_results/hdock_preds_{parental}/*.pdb")
        if DEBUG_MODE:
            paths_to_round1_hdock_predicted_poses = paths_to_round1_hdock_predicted_poses[0:3]
        main(
            paths_to_poses_list=paths_to_round1_hdock_predicted_poses,
            cfps_data_path=cfps_data_path,
            parental_seq=parental_seq,
            save_correlations_csv_path=save_correlations_csv_path,
            make_scatter_plots=True,
            save_plots_dir=f"{parental_results_dir}/plots/"
        )

        # tmux new -s pillbox
        # python3 run_pipeline.py 