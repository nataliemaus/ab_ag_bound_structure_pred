import pandas as pd 
import glob 
import numpy as np 
import os 
import copy 
import argparse 
from utils.get_mutations import get_mutations
from utils.refine_pose import refine_pose 
from utils.mutate_pose import mutate_residues 
from utils.score_pose import get_score_function
from utils.get_spearman_coeff import get_spearman_r
from utils.scatter_plot import scatter_plot
from constants import PARENTAL_ID_TO_AB_SEQ, HER2_CDR3_FIRST_INDEX, HER2_CDR3_LAST_INDEX, HER2_CDR_3


def main(
    hdock_pose_path,
    parental,
    affinity_data_path,
    parental_seq,
    results_dir,
    affinity_data_seq_col="HCDR3",
    affinity_data_label_col="-log(KD (M))",
    skip_refinement=False,
    mutations_only=True,
    save_pdbs=False, # new addition so I can send Ani intermediate PDBs to send through new model
):
    hdock_pose_number = int(hdock_pose_path.split("/")[-1].split("_")[-1].replace(".pdb", ""))
    # save_filename = hdock_pose_path.split("/")[-1].replace(".pdb", ".csv")
    save_filename = f"hdock_model_{hdock_pose_number}_results.csv"
    save_data_path = f"{results_dir}/{save_filename}"
    if save_pdbs:
        save_pdbs_dir = f"{results_dir}/hdock_model_{hdock_pose_number}_pdbs/"
        if not os.path.exists(save_pdbs_dir):
            os.mkdir(save_pdbs_dir) 

    df = pd.read_csv(affinity_data_path) # (422, 7)
    len_known_cdr3 = len(HER2_CDR_3)
    if mutations_only:
        # only keep new seqs that mutate original --> same length 
        len_new_cdr3s = [len(seq) for seq in df[affinity_data_seq_col].values]
        len_new_cdr3s = np.array(len_new_cdr3s)
        bool_array = len_new_cdr3s == len_known_cdr3
        df = df[bool_array] # (201,7)
    else:
        assert 0, "code not prepped to handle insertions or deletions"

    
    new_cdr3s = df[affinity_data_seq_col].values # .tolist()
    seqs_list = []
    for new_cdr3 in new_cdr3s:
        new_seq = copy.deepcopy(parental_seq)
        new_seq = [char for char in new_seq]
        new_cdr3 = [char for char in new_cdr3]
        new_seq[HER2_CDR3_FIRST_INDEX:HER2_CDR3_LAST_INDEX] = new_cdr3 
        new_seq = "".join(new_seq)
        seqs_list.append(new_seq)
    

    affinity_per_seq = df[affinity_data_label_col].values.astype(np.float32) # (201,)
    affinity_per_seq = affinity_per_seq.tolist()

    mutatant_positions_per_seq, mutatant_aas_per_seq = get_mutations(
        parental_seq=parental_seq,
        seqs_list=seqs_list,
    )

    energy_sfxn = get_score_function()

    binding_energies_per_seq = []
    mutated_pdb_paths_per_seq = []
    seq_number = 0
    for ix, mutatnt_positions_list in enumerate(mutatant_positions_per_seq):
        seq_number += 1
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
        # get binding energy
        binding_energy = energy_sfxn(refined_pose) 
        # save mutated and refined pose as pdb 
        if save_pdbs:
            save_pdb_filename = f"hdock{hdock_pose_number}_mutated_pose_{seq_number}.pdb"
            refined_pose.dump_pdb(f"{save_pdbs_dir}{save_pdb_filename}") 
            # TODO: save redfined pose 
            mutated_pdb_paths_per_seq.append(save_pdb_filename)
        
        # update and save data collected 
        binding_energies_per_seq.append(binding_energy)
        n_so_far = len(binding_energies_per_seq)
        be_array_so_far = np.array(binding_energies_per_seq)
        affinities_so_far = np.array(affinity_per_seq[0:n_so_far])
        if n_so_far > 2:
            spearman_r_so_far = get_spearman_r(affinities_so_far,be_array_so_far)
        else:
            spearman_r_so_far = -100
        # save to csv every loop so we have everything we've done so far if killed early 
        save_df = {}
        save_df["rosetta-energy"] = be_array_so_far
        save_df["known-affinity"] = affinities_so_far
        save_df["spearman_r"] = np.array([spearman_r_so_far]*n_so_far)
        if save_pdbs:
            save_df["mutated_pdb"] = np.array(mutated_pdb_paths_per_seq)
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
        help=" if True, do not do any refinement (refinement is the slowest step, skip to quick testing if env/setup works)",
        type=str2bool,
        default=False,
        required=False,
    ) 
    parser.add_argument(
        "--organize_data",
        help=" if True, run data organization rather than main",
        type=str2bool,
        default=False,
        required=False,
    ) 
    parser.add_argument(
        "--save_pdbs",
        help=" if True, run data organization rather than main",
        type=str2bool,
        default=False,
        required=False,
    ) 
    args = parser.parse_args() 

    parental = "HER2"
    if args.skip_refinement:
        results_dir = "her2/shzz_results_no_refinement"
    else:
        results_dir = "her2/shzz_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir) 
    

    if not args.organize_data: # run main with pose num 
        parental_seq = PARENTAL_ID_TO_AB_SEQ[parental]
        affinity_data_path = f"her2/Shanehsazzadeh-zero-shot-binders.csv"
        assert args.hdock_pose_num in np.arange(1,101).tolist()
        hdock_pose_path = f"her2/hdock_her2/model_{args.hdock_pose_num}.pdb"
        main(
            hdock_pose_path=hdock_pose_path,
            parental=parental,
            affinity_data_path=affinity_data_path,
            parental_seq=parental_seq,
            results_dir=results_dir,
            skip_refinement=args.skip_refinement,
            affinity_data_seq_col="HCDR3",
            affinity_data_label_col="-log(KD (M))",
            save_pdbs=args.save_pdbs,
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
    
    # NOTE: significant changes needed to run this for seqs with edits outside cdr3... 
    
    #   NOTE: can also run organize_data in the middle to see what it looks like before killing... 


    # NOTE: pose #3 is easily best and should get best correlation hopefully 

    # python3 shzz_data_run_pipeline_her2_parallel.py --organize_data False --hdock_pose_num 10
    # python3 shzz_data_run_pipeline_her2_parallel.py --organize_data True  
    
