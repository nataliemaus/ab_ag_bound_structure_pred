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

tmux attach -t ab1 - 20 

runai submit ab20 -v /home/nmaus/ab_ag_bound_structure_pred:/home/nmaus/ab_ag_bound_structure_pred --working-dir /home/nmaus/ab_ag_bound_structure_pred -i haydnj/torch:bayes-lqo -g 0 -e WANDB_API_KEY=fa9b0336bc46ee548faf75673c4f4ec5b461edb4 --interactive --attach

runai attach ab1

pip install Bio
pip3 install pyrosetta-installer 
python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'


# RUNNING 
# ab1-9 -->> poses 1-100 (redo)
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 1 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 2 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 3 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 4
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 5 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 6 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 7 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 8 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 9 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data False --hdock_pose_num 10 

# TODO: debug
# ab10-20 -->> poses 1-11 
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement False --organize_data False --hdock_pose_num 11
# all failed 
Traceback (most recent call last):
  File "/home/nmaus/ab_ag_bound_structure_pred/shzz_data_run_pipeline_her2_parallel.py", line 196, in <module>
    main(
  File "/home/nmaus/ab_ag_bound_structure_pred/shzz_data_run_pipeline_her2_parallel.py", line 79, in main
    refined_pose = refine_pose(
  File "/home/nmaus/ab_ag_bound_structure_pred/utils/refine_pose.py", line 298, in refine_pose
    os.remove(path_constraint_cst_file) 
FileNotFoundError: [Errno 2] No such file or directory: 'temp_constraints.cst'


python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement True --organize_data True  
python3 shzz_data_run_pipeline_her2_parallel.py --skip_refinement False --organize_data True  




# ------------------------
# DONE: 
python3 run_pipeline_her2_parallel.py --organize_data False --hdock_pose_num 2

python3 run_pipeline_her2_parallel.py --organize_data True 
