import pandas as pd 
import glob 
import numpy as np 
import os 
import copy 
import argparse 
from utils.mutate_utils import get_mutations_convert_seq1_to_seq2
from utils.refine_pose import refine_pose 
from utils.mutate_pose import mutate_residues 
from utils.score_pose import get_score_function
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot
from constants import PARENTAL_ID_TO_AB_SEQ, HER2_CDR3_FIRST_INDEX, HER2_CDR3_LAST_INDEX

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
    hdock_pose_path,
    parental,
    affinity_data_path,
    parental_seq,
    results_dir,
    skip_refinement=True,
):
    assert 0, "make sure to switch cdr3 in constants back to length 10 one for this dataset!"
    save_filename = hdock_pose_path.split("/")[-1].replace(".pdb", ".csv")
    save_data_path = f"{results_dir}/{save_filename}"

    df = pd.read_csv(affinity_data_path)

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

    binding_energies_per_seq = []
    for ix, mutatnt_positions_list in enumerate(mutatant_positions_per_seq):
        mutated_pose = mutate_residues(
            path_to_pdb=hdock_pose_path,
            parental=parental,
            mutant_positions=mutatnt_positions_list,
            mutant_aas=mutatant_aas_per_seq[ix],
            save_mutated_pdb=False,
        )
        if skip_refinement:
            refined_pose = mutated_pose
        else:
            refined_pose = refine_pose(
                parental=parental,
                pose=mutated_pose,
            ) 
        binding_energy = energy_sfxn(refined_pose) 
        binding_energies_per_seq.append(binding_energy)

        n_so_far = len(binding_energies_per_seq)
        be_array_so_far = np.array(binding_energies_per_seq)
        affinities_so_far = np.array(affinity_per_seq[0:n_so_far])
        if n_so_far > 2:
            spearman_r_so_far = get_spearman_r(affinities_so_far,be_array_so_far)
        else:
            spearman_r_so_far = -100

        save_df = {}
        save_df["energy"] = be_array_so_far
        save_df["affinity"] = affinities_so_far
        save_df["spearman_r"] = np.array([spearman_r_so_far]*n_so_far)
        save_df = pd.DataFrame.from_dict(save_df)
        save_df.to_csv(save_data_path, index=None)


def organize_data(
    hdock_pose_paths_list, 
    results_dir, 
    save_correlations_csv_path,
    make_scatter_plots=True,
):
    # for each pose in list, compute kd-energy spearman correlation coeff 
    save_plots_dir = f"{results_dir}/plots/"
    correlation_coeffs = []
    run_pose_paths = [] # paths to poses we've actually run main on to collect data 
    for hdock_pose_path in hdock_pose_paths_list:
        pose_data_path = hdock_pose_path.split("/")[-1].replace(".pdb", ".csv")
        pose_data_path = f"{results_dir}/{pose_data_path}"
        if os.path.exists(pose_data_path):
            pose_df = pd.read_csv(pose_data_path)
            binding_energies_per_seq = pose_df["energy"].values 
            affinity_per_seq = pose_df["affinity"].values 
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
        

    save_df = {}
    save_df["kd_energy_spearman_r"] = np.array(correlation_coeffs)
    save_df["pose_path"] = np.array(run_pose_paths)
    save_df = pd.DataFrame.from_dict(save_df)
    save_df = save_df.sort_values(by=["kd_energy_spearman_r"], ascending=False)
    save_df.to_csv(save_correlations_csv_path, index=None)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hdock_pose_num",
        help="Which hdock model to run.",
        type=int,
        default=1,
        required=False,
    )
    parser.add_argument(
        "--skip_refinement",
        help=" if True, do not do any refinement",
        type=str2bool,
        default=True,
        required=False,
    ) 
    parser.add_argument(
        "--organize_data",
        help=" if True, run data organization rather than main",
        type=str2bool,
        default=False,
        required=False,
    ) 
    args = parser.parse_args() 

    parental = "HER2"
    if args.skip_refinement:
        results_dir = "her2/results_no_refinement"
    else:
        results_dir = "her2/results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir) 
    
    
    if not args.organize_data: # run main with pose num 
        parental_seq = PARENTAL_ID_TO_AB_SEQ[parental]
        affinity_data_path = f"her2/affinity_data_cdr3.csv"
        assert args.hdock_pose_num in np.arange(1,101).tolist()
        hdock_pose_path = f"her2/hdock_her2/model_{args.hdock_pose_num}.pdb"
        main(
            hdock_pose_path=hdock_pose_path,
            parental=parental,
            affinity_data_path=affinity_data_path,
            parental_seq=parental_seq,
            results_dir=results_dir,
            skip_refinement=args.skip_refinement,
        )
    else: # run organize data 
        # get all hdock pdb paths 
        paths_to_round1_hdock_predicted_poses = glob.glob(f"her2/hdock_her2/model_*.pdb")
        # remove .clean.pdb versions 
        temp = []
        for hdock_pdb_file in paths_to_round1_hdock_predicted_poses:
            if not ("clean" in hdock_pdb_file):
                temp.append(hdock_pdb_file)
        paths_to_round1_hdock_predicted_poses = temp 
        # csv to save poses ordered by spearman_r
        save_correlations_csv_path = f"{results_dir}/{parental}_hdock_round1_poses_correlation_results.csv"
        make_scatter_plots = True 
        organize_data(
            hdock_pose_paths_list=paths_to_round1_hdock_predicted_poses, 
            results_dir=results_dir, 
            save_correlations_csv_path=save_correlations_csv_path,
            make_scatter_plots=make_scatter_plots,
        )
    
    # TODO: run for a long time for each of the 100 poses, kill when you feel like you have enought data 
    # then run organize_data to see what's up 
    #   NOTE: can also run organize_data in the middle to see what it looks like before killing... 


    # pose nums 1-10 running on Locust # tmux new -s ab1 , ... ab10 
    # NOTE: goal is to get 10-20 or so we don't need all 100, just has to include #3 ... 
    # python3 run_pipeline_her2_parallel.py --organize_data False --hdock_pose_num 10
    # python3 run_pipeline_her2_parallel.py --organize_data True  
    
