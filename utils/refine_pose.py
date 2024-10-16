# Some code copied from here: (igfold relax/refine)
# https://github.com/Graylab/IgFold/blob/main/igfold/refine/pyrosetta_ref.py
# 
# Constraint code grabbed from: 
# https://forum.rosettacommons.org/node/10105
# 
import os 
import pyrosetta 
import glob 
import pandas as pd 
import sys 
sys.path.append("../")
from utils.load_pose import load_pose 
from constants import PARENTAL_ID_TO_AG_SEQ, aa_1_to_3, aa_3_to_1, PARENTAL_ID_TO_CDR_INDEXES
from pyrosetta import rosetta
from rosetta.protocols.rigid import RigidBodyTransMover
import numpy as np 

init_string = "-mute all -ignore_zero_occupancy false -detect_disulf true -detect_disulf_tolerance 1.5 -check_cdr_chainbreaks false"
pyrosetta.init(init_string, silent=True)

def get_min_mover(
    max_iter: int = 1000,
    sf_name: str = "ref2015_cst",
    coord_cst_weight: float = 1, #
    dih_cst_weight: float = 1,
    constraint_applied: bool = True,
    constraint_weight: float = 1.0,
) -> pyrosetta.rosetta.protocols.moves.Mover:
    """
    Create full-atom minimization mover
    """

    sf = pyrosetta.create_score_function(sf_name)
    sf.set_weight(
        pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded,
        1,
    )
    sf.set_weight(
        pyrosetta.rosetta.core.scoring.ScoreType.pro_close,
        0,
    )
    sf.set_weight(
        pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint,
        coord_cst_weight,
    )
    sf.set_weight(
        pyrosetta.rosetta.core.scoring.ScoreType.dihedral_constraint,
        dih_cst_weight,
    )

    if constraint_applied: 
        score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        constraint = score_manager.score_type_from_name('atom_pair_constraint')
        sf.set_weight(constraint, constraint_weight)

    mmap = pyrosetta.rosetta.core.kinematics.MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(False)
    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover(
        mmap,
        sf,
        'lbfgs_armijo_nonmonotone',
        0.0001,
        True,
    )
    min_mover.max_iter(max_iter)
    min_mover.cartesian(True)

    return min_mover



def get_repack_mover():
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())

    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
    )
    packer.task_factory(tf)

    return packer


def write_list_of_lines_to_file(filepath, list_of_lines):
    f = open(filepath, "a")
    for line_w in list_of_lines:
        f.write(line_w )
    f.close()


def write_constraint_cst_file_from_pose(
    pose, 
    parental,
    save_path_constraint_cst_file, 
    harmonic_std=0.0001,
):
    # if this file already exists, delete and re-write it 
    if os.path.exists(save_path_constraint_cst_file):
        os.remove(save_path_constraint_cst_file)

    ag_seq = PARENTAL_ID_TO_AG_SEQ[parental]

    # ab and ag 
    # 2 chains for pillbox, 3 chains for her2 (her2 has 2 ab chains)
    chain_nums = {
        "ab":[],
    }
    for chain_num in range(1, pose.num_chains()+1):
        first_aa_in_chain = pose.residue(pose.chain_begin(chain_num)).name()
        if first_aa_in_chain.startswith(aa_1_to_3[ag_seq[0]]): # if chian is ag 
            chain_nums["ag"] = chain_num  
        else:
            chain_nums["ab"].append(chain_num) 

    if not ("ag" in chain_nums):
        import pdb 
        pdb.set_trace()
        # I think this happened bc I accidently tried to run on lig.pdb file!
        # TODO: investiage why this happened... 

        # (Pdb) chain_nums
        # {'ab': [1, 2]}
        # (Pdb) pose.num_chains()
        # 2
        # wtfff 

    # Used to get aa seqs for each chain for her2... 
    # aa_seqs = []
    # for chain_num in chain_nums["ab"]:
    #     aa_seq = ""
    #     for res_num in range(pose.chain_begin(chain_num), pose.chain_end(chain_num)+1):
    #         res1 = pose.residue(res_num).name()
    #         res1 = res1.split(":")[0]
    #         aa_1 = aa_3_to_1[res1]
    #         aa_seq += aa_1 
    #     aa_seqs.append(aa_seq)



    # first_aa_chain_1 = pose.residue(pose.chain_begin(1)).name()
    # first_aa_chain_2 = pose.residue(pose.chain_begin(2)).name()

    
    # if first_aa_chain_1.startswith(aa_1_to_3[ag_seq[0]]):
    #     ag_chain_num = 1
    #     ab_chain_num = 2
    # else:
    #     ag_chain_num = 2
    #     ab_chain_num = 1
    #     assert first_aa_chain_2.startswith(aa_1_to_3[ag_seq[0]])

    

    cdr_chain_id, cdr_indexes = PARENTAL_ID_TO_CDR_INDEXES[parental]

    for ix, ab_chain_num in enumerate(chain_nums["ab"]):
        ab_res_num_start = pose.chain_begin(ab_chain_num)
        if ix == cdr_chain_id:
            ab_cdr_nums = [idx + ab_res_num_start for idx in cdr_indexes]
        else:
            ab_cdr_nums = []
        ab_res_num_end = pose.chain_end(ab_chain_num)
        constraint_lines = []

        # add constraints to non-cdr ab residues: 
        pairs_added = []
        for res_n1 in range(ab_res_num_start, ab_res_num_end+1):
            res1 = pose.residue(res_n1)
            for res_n2 in range(ab_res_num_start, ab_res_num_end+1):
                # if not cdr
                if (res_n1 not in ab_cdr_nums) and (res_n2 not in ab_cdr_nums):
                    # if not same residue 
                    if (res_n1 != res_n2):
                        # if not constraint already added 
                        if (res_n1, res_n2) not in pairs_added:
                            res2 = pose.residue(res_n2)
                            # just constraint between first pair of atoms 
                            atom_name1 = res1.atom_name(1).replace(" ", "") # remove white space 
                            atom_name2 = res2.atom_name(1).replace(" ", "") 
                            atom1_location = res1.atom(1).xyz() 
                            atom2_location = res2.atom(1).xyz() 
                            dst_pair = atom1_location.distance(atom2_location) 
                            c_line = f"AtomPair {atom_name1} {res_n1} {atom_name2} {res_n2} HARMONIC {dst_pair} {harmonic_std}\n"
                            constraint_lines.append(c_line)
                            pairs_added.append((res_n1, res_n2))
                            pairs_added.append((res_n2, res_n1))


    # Add constraints to all residues in ag (don't allow ag to move): 
    ag_res_num_start = pose.chain_begin(chain_nums["ag"])
    ag_res_num_end = pose.chain_end(chain_nums["ag"])
    pairs_added = []
    for res_n1 in range(ag_res_num_start, ag_res_num_end+1):
        res1 = pose.residue(res_n1)
        for res_n2 in range(ag_res_num_start, ag_res_num_end+1):
            if (res_n1 != res_n2):
                if (res_n1, res_n2) not in pairs_added:
                    res2 = pose.residue(res_n2)
                    # just constrain between first pair of atoms 
                    atom_name1 = res1.atom_name(1).replace(" ", "") 
                    atom_name2 = res2.atom_name(1).replace(" ", "") 
                    atom1_location = res1.atom(1).xyz() 
                    atom2_location = res2.atom(1).xyz() 
                    dst_pair = atom1_location.distance(atom2_location) 
                    c_line = f"AtomPair {atom_name1} {res_n1} {atom_name2} {res_n2} HARMONIC {dst_pair} {harmonic_std}\n"
                    constraint_lines.append(c_line)
                    pairs_added.append((res_n1, res_n2))
                    pairs_added.append((res_n2, res_n1))

    # write constraint lines to cst file 
    write_list_of_lines_to_file(
        filepath=save_path_constraint_cst_file, 
        list_of_lines=constraint_lines,
    )


def apply_constraint(pose, constraints_cst_file_path):
    constraints = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(constraints_cst_file_path)
    constraints.add_constraints(True)
    constraints.apply(pose)


def refine(
    input_pdb_file,
    parental,
    out_pdb_file=None,
    minimization_iter=100,
    reuse_existing_pdbs=True,
    refine_cdrs_only=True,
    harmonic_std=0.0001, 
    delete_input_pdb_after_refine=True,
):  
    if out_pdb_file is None:
        out_pdb_file = input_pdb_file.replace(".pdb", f"_refined.pdb") 

    if reuse_existing_pdbs and os.path.exists(out_pdb_file):
        # if this path already exists, this refinement has arleady been done and saved, skip 
        return out_pdb_file

    pose = load_pose(input_pdb_file,) 

    if refine_cdrs_only:
        path_constraint_cst_file = input_pdb_file.replace(".pdb", f"_constraints.cst")
        write_constraint_cst_file_from_pose(
            pose=pose, 
            parental=parental,
            save_path_constraint_cst_file=path_constraint_cst_file,
            harmonic_std=harmonic_std,
        )
        apply_constraint(
            pose=pose, 
            constraints_cst_file_path=path_constraint_cst_file,
        )
        # remove cst file after constraint applied to save space 
        os.remove(path_constraint_cst_file) 

    min_mover = get_min_mover(
        max_iter=minimization_iter,
        coord_cst_weight=1,
        dih_cst_weight=0,
        constraint_applied=refine_cdrs_only,
    )
    min_mover.apply(pose)

    packer = get_repack_mover()
    packer.apply(pose)

    min_mover.apply(pose)

    pose.dump_pdb(out_pdb_file)

    # used to save space, only keep refined version 
    if delete_input_pdb_after_refine:
        os.remove(input_pdb_file)
    # always remove clean pdb file to save space 
    clean_pdb_file = input_pdb_file.replace(".pdb", ".clean.pdb")
    os.remove(clean_pdb_file)

    return out_pdb_file



def refine_pose(
    pose,
    parental,
    minimization_iter=100,
    refine_cdrs_only=True,
    harmonic_std=0.0001, 
):  
    if refine_cdrs_only:
        path_constraint_cst_file = "temp_constraints.cst"
        write_constraint_cst_file_from_pose(
            pose=pose, 
            parental=parental,
            save_path_constraint_cst_file=path_constraint_cst_file,
            harmonic_std=harmonic_std,
        )
        apply_constraint(
            pose=pose, 
            constraints_cst_file_path=path_constraint_cst_file,
        )
        # remove cst file after constraint applied to save space 
        os.remove(path_constraint_cst_file) 

    min_mover = get_min_mover(
        max_iter=minimization_iter,
        coord_cst_weight=1,
        dih_cst_weight=0,
        constraint_applied=refine_cdrs_only,
    )
    min_mover.apply(pose)

    packer = get_repack_mover()
    packer.apply(pose)

    min_mover.apply(pose)

    return pose 

