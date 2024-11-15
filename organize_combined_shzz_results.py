import numpy as np 
import pandas as pd 
import os 
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot
import glob 

def kd_nM_to_neglog_kd_M(kd_nM):
    # Input: array of kds in raw nM ("KD (nM)" col in shzz data)
    # Output: array converged -log(KD (M)) ("-log(KD (M))" col in shzz data)
    kd_M = kd_nM*1e-9 
    neglog_kd_M = -np.log10(kd_M) 

    return neglog_kd_M


def organize_data(
    skip_refinement=True,
    take_neg_log=True,
):
    make_scatter_plots=True
    og_zero_shot_df = pd.read_csv("her2/Shanehsazzadeh-zero-shot-binders.csv")
    # HCDR3,KD (nM),-log(KD (M)),
    og_zero_shot_kds = og_zero_shot_df["KD (nM)"].values.astype(np.float32).tolist()
    og_zero_shot_neglog_kds = og_zero_shot_df["-log(KD (M))"].values.astype(np.float32).tolist()

    paths_to_round1_hdock_predicted_poses = glob.glob(f"her2/hdock_her2/model_*.pdb")
    # for each pose in list, compute kd-energy spearman correlation coeff 
    control_data_dir = "her2/control_shzz_results"
    zero_shot_data_dir = "her2/shzz_results"
    save_results_dir = "her2/zero_control_combined_shzz_results"
    if take_neg_log:
        save_results_dir = save_results_dir + "_neglogkd"
    if skip_refinement:
        control_data_dir = control_data_dir + "_no_refinement"
        zero_shot_data_dir = zero_shot_data_dir + "_no_refinement"
        save_results_dir = save_results_dir + "_no_refinement"
    if not os.path.exists(save_results_dir):
        os.mkdir(save_results_dir)
    save_correlations_csv_path = f"{save_results_dir}/HER2_hdock_round1_poses_correlation_results.csv"
    save_plots_dir = f"{save_results_dir}/plots/"
    correlation_coeffs = []
    run_pose_paths = [] # paths to poses we've actually run main on to collect data 
    for hdock_pose_path in paths_to_round1_hdock_predicted_poses:
        pose_data_filename = hdock_pose_path.split("/")[-1].replace(".pdb", ".csv")

        control_data_path = f"{control_data_dir}/{pose_data_filename}"
        zero_shot_data_path = f"{zero_shot_data_dir}/{pose_data_filename}"
        # pose_data_path = f"{results_dir}/{pose_data_filename}"
        if os.path.exists(control_data_path) and os.path.exists(zero_shot_data_path):
            control_df = pd.read_csv(control_data_path)
            control_be = control_df["energy"].values 
            control_kds = control_df["affinity"].values 

            zero_df = pd.read_csv(zero_shot_data_path)
            zero_be = zero_df["energy"].values 
            zero_neglog_kds = zero_df["affinity"].values 
            zero_kds = []
            for neglog_kd in zero_neglog_kds:
                og_index = og_zero_shot_neglog_kds.index(neglog_kd) # will throw error if not in list 
                zero_kds.append(og_zero_shot_kds[og_index])
            zero_kds = np.array(zero_kds)

            affinity_per_seq = np.array(control_kds.tolist() + zero_kds.tolist() )
            if take_neg_log:
                affinity_per_seq  = kd_nM_to_neglog_kd_M(affinity_per_seq )
            binding_energies_per_seq = np.array(control_be.tolist() + zero_be.tolist() )
            # TODO: LATER? save combined df with marker col for control vs zero ?
            spearman_r = get_spearman_r(affinity_per_seq,binding_energies_per_seq)
            if make_scatter_plots:
                plot_filename = hdock_pose_path.split("/")[-1].replace(".pdb", ".png")
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
            run_pose_paths.append(hdock_pose_path)

    ascending = not take_neg_log
    save_df = {}
    save_df["kd_energy_spearman_r"] = np.array(correlation_coeffs)
    save_df["pose_path"] = np.array(run_pose_paths)
    save_df = pd.DataFrame.from_dict(save_df)
    save_df = save_df.sort_values(by=["kd_energy_spearman_r"], ascending=ascending)
    save_df.to_csv(save_correlations_csv_path, index=None)

if __name__ == "__main__":
    organize_data(
        skip_refinement=True,
        take_neg_log=True,
    )