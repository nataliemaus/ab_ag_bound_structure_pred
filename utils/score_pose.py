import pyrosetta 
pyrosetta.init()
from pyrosetta.teaching import get_fa_scorefxn


def get_score_function():
    scorefxn = get_fa_scorefxn()
    return scorefxn


def score_pose(pose):
    scorefxn = get_score_function()
    return scorefxn(pose)
