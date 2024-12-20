import os 
import pandas as pd 
import numpy as np 
import glob 
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot


affinity_keys = ["h3_mean"]# , "h1_mean","fluB_mean"] # ,"h1_label","h3_label","fluB_label"]

def organize_v2_all_affinity_data_and_seqs(
    results_dir, 
    make_scatter_plots=True,
    aff_percentile=None,
    aff_min_val=6.1,
):
    for affinity_key in affinity_keys:
        # for each pose in list, compute kd-energy spearman correlation coeff 
        affinity_key_results_dir = f"{results_dir}/{affinity_key}"
        if not os.path.exists(affinity_key_results_dir):
            os.mkdir(affinity_key_results_dir) 

        save_correlations_csv_path = f"{affinity_key_results_dir}/{affinity_key}_correlations.csv"
        save_plots_dir = f"{affinity_key_results_dir}/plots/"
        if not os.path.exists(save_plots_dir):
            os.mkdir(save_plots_dir) 
        correlation_coeffs = []
        hodck_pose_nums = [] # paths to poses we've actually run main on to collect data 
        for i in range(100):
            hdock_pose_num = i+1 
            pose_data_path = f"{results_dir}/model_{hdock_pose_num}.csv"
            pose_df = pd.read_csv(pose_data_path)
            binding_energies_per_seq = pose_df["energy"].values 
            affinity_per_seq = pose_df[affinity_key].values 
            bool_arr1 = np.logical_not(np.isnan(binding_energies_per_seq))
            bool_arr2 = np.logical_not(np.isnan(affinity_per_seq))
            bool_arr3 = np.isreal(binding_energies_per_seq)
            bool_arr4 = np.isreal(affinity_per_seq)
            bool_arr = np.logical_and(bool_arr1, bool_arr2)
            bool_arr = np.logical_and(bool_arr, bool_arr3)
            bool_arr = np.logical_and(bool_arr, bool_arr4)
            binding_energies_per_seq = binding_energies_per_seq[bool_arr]
            affinity_per_seq = affinity_per_seq[bool_arr]
            if aff_percentile is not None:
                min_val = np.percentile(affinity_per_seq, aff_percentile)
                bool_arr_high_aff = affinity_per_seq > min_val
                binding_energies_per_seq = binding_energies_per_seq[bool_arr_high_aff]
                affinity_per_seq = affinity_per_seq[bool_arr_high_aff]
            if aff_min_val is not None:
                bool_arr_high_aff = affinity_per_seq > aff_min_val
                binding_energies_per_seq = binding_energies_per_seq[bool_arr_high_aff]
                affinity_per_seq = affinity_per_seq[bool_arr_high_aff]


            spearman_r = get_spearman_r(affinity_per_seq,binding_energies_per_seq)
        
            if make_scatter_plots:
                plot_filename = f"model_{hdock_pose_num}.png"
                rounded_spearman_r = round(spearman_r, 3)
                scatter_plot( 
                    x=binding_energies_per_seq, 
                    y=affinity_per_seq, 
                    x_label="Rosetta Binding Energy", 
                    y_label=f"Affinity {affinity_key}", 
                    title=f"Affinity-Energy Correlation Spearman Coeff = {rounded_spearman_r}", 
                    path_to_save_plot=f"{save_plots_dir}{plot_filename}", 
                )
            correlation_coeffs.append(spearman_r)
            hodck_pose_nums.append(hdock_pose_num)

        save_df = {}
        save_df["kd_energy_spearman_r"] = np.array(correlation_coeffs)
        save_df["hodck_pose_num"] = np.array(hodck_pose_nums)
        save_df = pd.DataFrame.from_dict(save_df)
        save_df = save_df.sort_values(by=["kd_energy_spearman_r"], ascending=False)
        save_df.to_csv(save_correlations_csv_path, index=None)


if __name__ == "__main__":
    parental = "CR9114"
    results_dir = f"influenza/{parental}_all_seqs_and_affinity_data_no_refinement"
    organize_v2_all_affinity_data_and_seqs(
        results_dir=results_dir, 
    )