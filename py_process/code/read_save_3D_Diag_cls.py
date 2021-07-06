#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:10:10 2021

@author: tpeng
"""


import sys
import os
home_dir=os.path.expanduser("~")
package_path=home_dir+'/MITgcm_post/py_process/ppmod/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate

#%% 


import argparse
import shutil


parser_diag_cls = argparse.ArgumentParser(description="Load Diagnostics")
parser_diag_cls.add_argument('--dx', metavar="dx", type=int, default=500, help="horizontal resolution (x-dir)")
parser_diag_cls.add_argument('--dy', metavar="dy", type=int, default=500, help="horizontal resolution (y-dir)")
parser_diag_cls.add_argument('--Lx', metavar="Lx", type=int, default=480e3, help="horizontal resolution (x-dir)")
parser_diag_cls.add_argument('--Ly', metavar="Ly", type=int, default=960e3, help="horizontal resolution (y-dir)")    
parser_diag_cls.add_argument('--ind_z', metavar="ind_z", default='1,2,3,4,5,20,60', help="vertical indices recorded (see data.diagnostics)")    
parser_diag_cls.add_argument('--fhead', metavar="fhead", help="Filename's head string of a field recorded (see data.diagnostics)")    
parser_diag_cls.add_argument('--path_scratch',metavar="path_scratch", help="where raw results are stored.")    
parser_diag_cls.add_argument('--path_results', metavar="path_results", help="where extracted results are saved.")    
parser_diag_cls.add_argument('--groupname', metavar="--groupname", help="Case name (casename)")
parser_diag_cls.add_argument('--casename', metavar="casename", help="Case name (casename)")
parser_diag_cls.add_argument('--start_ind', type=int, metavar="start_ind", help="Start index (start_ind)")
parser_diag_cls.add_argument('--end_ind', type=int, metavar="end_ind", help="Start index (end_ind)")


parser_diag_cls.add_argument('--tape_days', metavar="tape_days", type=int, help="how many days in a tape")
parser_diag_cls.add_argument('--Fs', metavar="Fs", type=int, help="how many points sampled per day (86400 s)")
parser_diag_cls.add_argument('--dt_model', metavar="dt_model", type=int, help="how many points sampled per day (86400 s)")
parser_diag_cls.add_argument('--z_target', metavar="z_target", help="which layer to extract (can be longer than 1")
parser_diag_cls.add_argument('--name_fields', metavar="name_fields", help="which layer to extract (can be longer than 1")

args_cls = parser_diag_cls.parse_args()

dict_load_diag={
    'dx' : int(args_cls.dx),
    'dy' : int(args_cls.dy),
    'Lx' : int(args_cls.Lx), 
    'Ly' : int(args_cls.Ly),
    'ind_z' : np.array((args_cls.ind_z).split(','),dtype=int),
    'fhead' : args_cls.fhead,
    'path_scratch' : args_cls.path_scratch,
    'path_results' : args_cls.path_results,
    'groupname' : args_cls.groupname,
    'casename' : args_cls.casename,
    'start_ind' : args_cls.start_ind,
    'end_ind' : args_cls.end_ind
    }
print('Parameters read from at runtime')
print(dict_load_diag)


dict_inside_cls={
    'tape_days' : args_cls.tape_days, # tape every xx days
    'Fs' : args_cls.Fs, # samples a day
    'dt_model' : args_cls.dt_model, # [sec] timestep in the simulations
    'z_target' : np.array(args_cls.z_target.split(','),dtype=int), # [1] for surface
    'name_fields' : args_cls.name_fields.split(',') # name used as the head of saved fields (tensor)
    }
print(dict_inside_cls)



# dict_load_diag={
#     'dx' : 500.0,
#     'dy' : 500.0,
#     'Lx' : 480e3, 
#     'Ly' : 960e3,
#     'ind_z' : [1,2,3,4,5,20,60],
#     'fhead' : 'Diag_snaps_UV',
#     'path_scratch' : home_dir+'/scratch/MITgcm_cases/',
#     'path_results' : home_dir+'/postproc/results/',
#     'groupname' : '',
#     'casename' : 'rel_uvar_9',
#     'start_ind' : 700000,
#     'end_ind' : 715000
#     }
# dict_inside_cls={
#     'tape_days' : 14, # tape every xx days
#     'Fs' : 8, # samples a day
#     'dt_model' : 75, # [sec] timestep in the simulations
#     'z_target' : [1], # [1] for surface
#     'name_fields' : ['Uvel', 'Vvel'] # name used as the head of saved fields (tensor)
#     }

#%%
def load_diag_files(Diag_obj,tape_days,Fs,dt_model,z_target,name_fields):
    Diag_obj.assign_path()
    Diag_obj.get_file_ind()
    Diag_obj.tape_file(tape_days,Fs,dt_model)
    Diag_obj.cut_xylayers([1])
    Diag_obj.load_diag_files(name_fields)
    return
#%%

Diag_UV_7layers=cls_3D_Diag(**dict_load_diag)
load_diag_files(Diag_UV_7layers,**dict_inside_cls)

#%% save object
for k in arange(len(dict_inside_cls['name_fields'])):
    path_save_cls=Diag_UV_7layers.path_results+'/saved_obj/'+Diag_UV_7layers.groupname+'/'+Diag_UV_7layers.casename+'/'
    if not os.path.exists(path_save_cls):
        os.mkdir(path_save_cls)
    file_save_cls = open(path_save_cls+dict_inside_cls['name_fields'][k]+'.obj', 'wb') 
    pickle.dump(Diag_UV_7layers,file_save_cls)
