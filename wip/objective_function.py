import numpy as np 
from utils.mutate_utils import get_mutations_convert_seq1_to_seq2
from utils.mutate_pose import mutate_residues
from utils.refine_pose import refine_pose  
from utils.score_pose import get_score_function 
from constants import PARENTAL_ID_TO_AA_SEQ


class AffinityStructureObjective():
    ''' Objective function for optimizing affinity and variance in binding energy across possible poses 
        Goal: acheive high affinity and high variance in binding energy among hdock poses 
    ''' 
    def __init__( 
        self,
        hdock_poses_list, # list of paths to hdock poses 
        parental_seq,
        affinity_weight=1.0,
        n_mutations_threshold=5,
        refinement_minimization_iter=100,
        variance_normalizing_constant=3_000,
    ):
        self.hdock_poses_list = hdock_poses_list 
        self.parental_seq = parental_seq 
        # self.bighat_affinity_oracle_wrapper = TODO 
        self.affinity_weight = affinity_weight
        self.n_mutations_threshold = n_mutations_threshold
        self.refinement_minimization_iter = refinement_minimization_iter
        self.energy_sfxn = get_score_function()
        self.variance_normalizing_constant = variance_normalizing_constant


    def bighat_affinity_oracle_wrapper(self, x_list):
        # TODO: replace with actual bighat affinity oracle wrapper, 
        # for now this is ued for testing 
        kds = []
        for x in x_list:
            kds.append(10.0)
        return kds 

    def query_black_box(self, x_list):
        assert len(x_list) >= 2
        affinities = self.bighat_affinity_oracle_wrapper(x_list)
        rewards = []
        for i, x_i in enumerate(x_list):
            binding_energy_variance = self.get_energy_variance(x_i)
            reward = (affinities[i]*self.affinity_weight) + (binding_energy_variance//self.variance_normalizing_constant)
            rewards.append(reward)

        return rewards 
    
    def get_energy_variance(self, x):
        # 1. get list of mutations done to create seq x 
        positions_list, aas_list, n_mutations = get_mutations_convert_seq1_to_seq2(
            seq1=self.parental_seq,
            seq2=x,
        )
        if n_mutations > self.n_mutations_threshold: 
            return -100

        # 1. get binding energy variance 
        energies = []
        for pose_path_i in self.hdock_poses_list:
            mutated_pose = mutate_residues(
                path_to_pdb=pose_path_i,
                mutant_positions=positions_list, 
                mutant_aas=aas_list, 
                save_mutated_pdb=False,
            )
            refined_pose = refine_pose(
                pose=mutated_pose,
                minimization_iter=self.refinement_minimization_iter,
            ) 
            binding_energy = self.energy_sfxn(refined_pose) 
            print(binding_energy)
            energies.append(binding_energy) 
        energy_var = np.var(np.array(energies)) 

        return energy_var


# TODO: add thermostability constraint :) 


if __name__ == "__main__":
    hdock_model_nums_list=[
        7, 51, 58, 82, 93, 26, 38, 100, 59, 2, 36, 84, 11, 1, 79, 96, 14, 15, 88, 89
    ]
    hdock_model_nums_list = hdock_model_nums_list[0:4] # smaller list of quick testing 
    hdock_poses_list = [
        f"hdock_preds_A38781/model_{hdock_num}.pdb" for hdock_num in hdock_model_nums_list
    ]
    obj = AffinityStructureObjective(
        hdock_poses_list=hdock_poses_list,
        parental_seq = PARENTAL_ID_TO_AA_SEQ["A38781"], 
        refinement_minimization_iter=10, # smaller n iters for quick testing 
    )
    seqs = [
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNDLKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASTNTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKWRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVLDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRKAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLCNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDKAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQCDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAARGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTIYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTIIRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYCGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGCYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYFWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNEKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYVVPVKGRFTISRDNANKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYDVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
        "QVQLVESGGGLVQAGGSLRLSCAASGSVSDTNLMGWYRQAPGKQRDLVASITTGGLRNYIVPVKGRFTISRDNAKKTSYLQMNALKPDDTAVYYCWSRDIFRGNEYWGQGTQVTVSS",
    ]
    rewards = obj.query_black_box(seqs)
    print(rewards)
