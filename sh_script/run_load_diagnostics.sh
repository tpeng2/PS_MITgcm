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
end_ind=718446

tape_days=7
Fs=8
dt_model=75
z_target=1

# ===========
# file #1
# ===========

name_fields='Uvel','Vvel'

python /home/tpeng/MITgcm_post/py_process/code/read_save_3D_Diag_cls.py --dx=$dx --dy=$dy\
      --Lx=$Lx --Ly=$Ly --ind_z=$ind_z --fhead=$fhead --path_scratch=$path_scratch\
      --path_results=$path_results --groupname=$groupname  --casename=$casename\
      --start_ind=$start_ind --end_ind=$end_ind --tape_days=$tape_days --Fs=$Fs\
      --dt_model=$dt_model --z_target=$z_target --name_fields=$name_fields
# =============
# file #2
# =============
fhead='Diag_snaps_WT'
name_fields='Wvel' #,'THETA'
python /home/tpeng/MITgcm_post/py_process/code/read_save_3D_Diag_cls.py --dx=$dx --dy=$dy\
      --Lx=$Lx --Ly=$Ly --ind_z=$ind_z --fhead=$fhead --path_scratch=$path_scratch\
      --path_results=$path_results --groupname=$groupname  --casename=$casename\
      --start_ind=$start_ind --end_ind=$end_ind --tape_days=$tape_days --Fs=$Fs\
      --dt_model=$dt_model --z_target=$z_target --name_fields=$name_fields

# =============
# file #3
# =============
fhead='PHIHYD'
name_fields='P'
python /home/tpeng/MITgcm_post/py_process/code/read_save_3D_Diag_cls.py --dx=$dx --dy=$dy\
      --Lx=$Lx --Ly=$Ly --ind_z=$ind_z --fhead=$fhead --path_scratch=$path_scratch\
      --path_results=$path_results --groupname=$groupname  --casename=$casename\
      --start_ind=$start_ind --end_ind=$end_ind --tape_days=$tape_days --Fs=$Fs\
      --dt_model=$dt_model --z_target=$z_target --name_fields=$name_fields


