# https://pypi.org/project/pyrosetta-installer/ 
pip3 install pyrosetta-installer 
python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

XXX 
pip install pyrosetta-installer 
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

# Above worked on allegro, but not on my laptop --> will need to run everything on allegro. 

#  git repo: 
https://github.com/nataliemaus/ab_ag_bound_structure_pred

# On Locust: 

tmux attach -t ab1 - 34

runai submit ab20 -v /home/nmaus/ab_ag_bound_structure_pred:/home/nmaus/ab_ag_bound_structure_pred --working-dir /home/nmaus/ab_ag_bound_structure_pred -i haydnj/torch:bayes-lqo -g 0 -e WANDB_API_KEY=fa9b0336bc46ee548faf75673c4f4ec5b461edb4 --interactive --attach

runai attach ab1

pip install Bio
pip3 install pyrosetta-installer 
python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

# ------------ SHZZ ZERO-SHOT -------------------

# DONE + DATA DOWNLOADED: poses 1-100  no refinement
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 100 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data True  

# DONE + DATA DOWNLOADED: poses 1-100 w/ refinement
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement False --organize_data False --hdock_pose_num 100
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement False --organize_data True 

# ------------ SHZZ Controls (more data) -------------------

# RUNNING: SHZZ poses 1-100 NO refinment (DONE, DOWNLOADED)
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 1
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement True --organize_data True 


# RUNNING: SHZZ poses 1-100 WITH refinment 
# ab1-32 x 3 runs each --> poses 1-100 running (queue order: 1-33, 34-66, 67-100)
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement False --organize_data False --hdock_pose_num 1
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement False --organize_data False --hdock_pose_num 33
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement False --organize_data False --hdock_pose_num 67 (accidently started on 66 --> repeat --> aabo34 running 98,99 )

# Then do: 
python3 control_shzz_run_pipeline_her2_parallel.py --skip_refinement False --organize_data True 



# ------------------------ RUNNING INFLUENZA -------------------------------


# Running: poses 1-100 without refinment first 
# V1 looked bad, now trying V2 where we ignore any mutations in first or last 6 aa's (V2), if this is bad maybe there's something else wrong going on? 
# ab1-20 (5 per gpu)
# Kicked off Tuesday Dec 3rd 7:25am 
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data False --skip_refinement True --hdock_pose_num 20
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data False --skip_refinement True --hdock_pose_num 40
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data False --skip_refinement True --hdock_pose_num 60
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data False --skip_refinement True --hdock_pose_num 80
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data False --skip_refinement True --hdock_pose_num 100


# Then do: 
python3 inf_run_pipeline_parallel.py --num_affinity_data 1000 --organize_data True --skip_refinement True




# Version with all seqs and save all affinity data:
#   running ab1-20 x 5 poses each 
python3 inf_run_pipeline_parallel.py --organize_data False --use_all_affinity_data_and_seqs True --skip_refinement True --hdock_pose_num 20
python3 inf_run_pipeline_parallel.py --organize_data False --use_all_affinity_data_and_seqs True --skip_refinement True --hdock_pose_num 40 
python3 inf_run_pipeline_parallel.py --organize_data False --use_all_affinity_data_and_seqs True --skip_refinement True --hdock_pose_num 60 
python3 inf_run_pipeline_parallel.py --organize_data False --use_all_affinity_data_and_seqs True --skip_refinement True --hdock_pose_num 80 
python3 inf_run_pipeline_parallel.py --organize_data False --use_all_affinity_data_and_seqs True --skip_refinement True --hdock_pose_num 100
 
# TODO NEXT
# 1. never version of organize data to compute spearman_r's and plots 
#       using different combos of affinity cols and ranges of values 
#       can we get good correlations for any? (i.e. only >6.0, 6.5)
#       can do locally! 