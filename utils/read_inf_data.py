import pandas as pd 
import numpy as np 


def load_influenza_affinity_data(
    inf_dir="../influenza",
    ab_parental="CR9114",
    ag_h_num=3,
    N_load=600,
    get_all_affinity_labels_and_seqs=False,
): 
    assert ab_parental == "CR9114", "hdock poses only computed for CR9114, not yet for CR6261"
    all_affinity_labels = ['h1_mean', 'h3_mean', 'fluB_mean', 'h1_label', 'h3_label', 'fluB_label']
    data_dir = f"{inf_dir}/zenodo_affinity_data/phillips_et_al_processed/{ab_parental}"
    df = pd.read_csv(f"{data_dir}/all_combined_data.csv") #  (65536,) 
    if get_all_affinity_labels_and_seqs:
        print(f"N_load={N_load} is irrelevant since get_all_affinity_labels_and_seqs=True")
        print(f"Loading all N={df.shape[0]} seqs and labels ...")
    else:
        assert ag_h_num == 3, "hdock poses only computed for CR9114 + H3 Influenza AG so far"

    # remove mutations (seqs not equal to length 121)
    lengths = np.array([len(seq) for seq in df['H_seq'].values]) # (600,)
    len_bool_arr = lengths == 121 
    df = df[len_bool_arr] # (65536,)
    seqs = df['H_seq'].values #  (65536,) 

    if get_all_affinity_labels_and_seqs: 
        affinity_labels_dict = {}
        for affinity_label in all_affinity_labels:
            affinity_labels_dict[affinity_label] = df[affinity_label].values
        return seqs, affinity_labels_dict

    
    # Otherwise use f"h{ag_h_num}_mean" labels and get roughly equal number per affinity bin... 
    affinities = df[f"h{ag_h_num}_mean"].values # (65536,)
    # remove nan affinities 
    bool_arr = np.logical_not(np.isnan(affinities))
    seqs = seqs[bool_arr] # (65535,)
    affinities = affinities[bool_arr] # (65535,)
    #   None removed, all have same length of 121 
    # Get a balanced dataset with good number in each affinity range 
    # NOTE: this is catered to the CR9114 + ag_h_num=3 dataset, re-examine bins for others
    bins = {}
    bins[1] = affinities > 8.5 # (8.5-9] 121 total seqs 
    bins[2] = np.logical_and(affinities > 8.0, affinities <= 8.5) # (8.0-8.5] 264 total seqs 
    bins[3] = np.logical_and(affinities > 7.0, affinities <= 8.0) # (7.0-8.0] 2881 total seqs 
    bins[4] = np.logical_and(affinities > 6.0, affinities <= 7.0) # (6.0-7.0] 3908 total seqs 
    bins[5] = affinities == 6.0 # (6.0, non-binders) 58361 total seqs 
    n_bins = len(bins.keys()) # 5
    n_seqs_per_bin = N_load//n_bins # 120
    # add an equal number of indexes from each bin 
    good_indexes = []
    for bin_num in bins.keys():
        good_indexes = good_indexes + np.where(bins[bin_num])[0][0:n_seqs_per_bin].tolist()
    affinities = affinities[good_indexes] # (600,)
    seqs = seqs[good_indexes] # (600,)

    
    return seqs, affinities 



    if False: # initial testing of data 
        # (Pdb) df.keys() 
        # Index(['binary_id', 'str_aa_id', 'edit_distance', 'H_seq', 'h1_mean',
        #    'h3_mean', 'fluB_mean', 'h1_label', 'h3_label', 'fluB_label'],
        #   dtype='object')
    
        affinities = df[f"h3_mean"].values # 
        affinities = affinities[np.logical_not(np.isnan(affinities))] # 65535
        # (Pdb) affinities.min()
        # 6.0
        # (Pdb) affinities.max()
        # 9.041930737842575
        # (Pdb) affinities.mean()
        # 6.106048829827008
        n_greater_than_6 = (affinities > 6.0).sum() # 7174
        n_greater_than_7 = (affinities > 7.0).sum() # 3266
        n_greater_than_8 = (affinities > 8.0).sum() # 385 
        n_greater_than_8p5 = (affinities > 8.5).sum() # 121
        n_greater_than_8p6 = (affinities > 8.6).sum() # 83
        n_greater_than_8p7 = (affinities > 8.7).sum() # 45 
        n_greater_than_8p8 = (affinities > 8.8).sum() # 20
        n_greater_than_8p9 = (affinities > 8.9).sum() # 3
        n_greater_than_9 = (affinities > 9.0).sum() # 2
        n_geq_9 = (affinities >= 9.0).sum() # 2


        affinities2 = df['fluB_mean'].values # 
        affinities2 = affinities2[np.logical_not(np.isnan(affinities2))]
        # (Pdb) affinities2.mean(), affinities2.min(), affinities2.max()
        # (6.003074421673384, 6.0, 8.060834302617124)
        n_greater_than_6 = (affinities2 > 6.0).sum() # 198


        affinities3 = df['h3_label'].values #
        affinities3 = affinities3[np.logical_not(np.isnan(affinities3))]
        # (Pdb) np.unique(affinities3)  array([0, 1])
        n_binders = (affinities3 == 1).sum() # 7174
    
    
    


if __name__ == "__main__":
    seqs, affinities = load_influenza_affinity_data() 
    print(seqs.shape, affinities.shape, seqs[0], affinities.max(), affinities.min())