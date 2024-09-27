import os 
import sys 
sys.path.append("../")
from constants import PARENTAL_ID_TO_AA_SEQ
import pyrosetta 
pyrosetta.init()


def load_pose(path_to_pdb):
    clean_pdb_path = path_to_pdb.replace(".pdb", ".clean.pdb")
    if os.path.exists(clean_pdb_path):
        pose = pyrosetta.io.pose_from_pdb(clean_pdb_path)
    else:
        # first clean pdb to get pyrosetta desired format, creates .clean.pdb version 
        pyrosetta.toolbox.cleanATOM(path_to_pdb)
        pose = pyrosetta.io.pose_from_pdb(clean_pdb_path)
    return pose 
