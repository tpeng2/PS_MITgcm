#!/bin/bash

dx=500
dy=500
Lx=480000
Ly=960000
ind_z=1,2,3,4,5,20,60
fhead='Diag_snaps_UV'
path_scratch=$HOME/scratch/MITgcm_cases/
#path_scratch=/mnt/e/
path_results=$HOME/postproc/results/
groupname=''
casename="rel_uvar_9"
start_ind=700000
end_ind=734421

tape_days=7
Fs=8
dt_model=75
z_target=1

# =============
# plotting
# =============

python $HOME/MITgcm_post/py_process/code/load_3D_Diag_plot_kw.py --dx=$dx --dy=$dy\
	--Lx=$Lx --Ly=$Ly --path_scratch=$path_scratch --path_results=$path_results\
	--groupname=$groupname  --casename=$casename --start_ind=$start_ind --end_ind=$end_ind\
	--tape_days=$tape_days --Fs=$Fs --dt_model=$dt_model --z_target=$z_target --name_fields=$name_fields

