import pandas as pd 

def load_inf_affinity_data(
    inf_dir="../influenza",
    inf_ab_parental="CR9114",
): 
    assert inf_ab_parental == "CR9114", "hdock poses only computed for CR9114, not yet for CR6261"
    data_dir = f"{inf_dir}/zenodo_affinity_data/phillips_et_al_processed/{inf_ab_parental}"
    df = pd.read_csv(f"{data_dir}/all_combined_data.csv")
    import pdb 
    pdb.set_trace()
    # (Pdb) df.keys() 
    # Index(['binary_id', 'str_aa_id', 'edit_distance', 'H_seq', 'h1_mean',
    #    'h3_mean', 'fluB_mean', 'h1_label', 'h3_label', 'fluB_label'],
    #   dtype='object')


if __name__ == "__main__":
    load_inf_affinity_data() 