#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 11:37:24 2021
Read extracted files and process them
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
import PPFCN.fftfcn as fcn

from scipy import interpolate
import pickle

#%%

import argparse
import shutil


parser_diag_cls = argparse.ArgumentParser(description="Load Diagnostics")
parser_diag_cls.add_argument('--dx', metavar="dx", type=int, default=500, help="horizontal resolution (x-dir)")
parser_diag_cls.add_argument('--dy', metavar="dy", type=int, default=500, help="horizontal resolution (y-dir)")
parser_diag_cls.add_argument('--Lx', metavar="Lx", type=int, default=480e3, help="horizontal resolution (x-dir)")
parser_diag_cls.add_argument('--Ly', metavar="Ly", type=int, default=960e3, help="horizontal resolution (y-dir)")    
parser_diag_cls.add_argument('--path_scratch',metavar="path_scratch", default=home_dir+'/scratch/MITgcm_cases/', help="where raw results are stored.")    
parser_diag_cls.add_argument('--path_results', metavar="path_results", default=home_dir+'/postproc/results/',help="where extracted results are saved.")    
parser_diag_cls.add_argument('--groupname', metavar="--groupname", help="Case name (casename)")
parser_diag_cls.add_argument('--casename', metavar="casename", help="Case name (casename)")
parser_diag_cls.add_argument('--start_ind', type=int, metavar="start_ind", help="Start index (start_ind)")
parser_diag_cls.add_argument('--end_ind', type=int, metavar="end_ind", help="Start index (end_ind)")


parser_diag_cls.add_argument('--tape_days', metavar="tape_days", type=int, help="how many days in a tape")
parser_diag_cls.add_argument('--Fs', metavar="Fs", type=int, help="how many points sampled per day (86400 s)")
parser_diag_cls.add_argument('--dt_model', metavar="dt_model", default=75, type=int, help="how many points sampled per day (86400 s)")
parser_diag_cls.add_argument('--z_target', metavar="z_target", help="which layer to extract (can be longer than 1")
parser_diag_cls.add_argument('--name_fields', metavar="name_fields", help="which layer to extract (can be longer than 1")


def plot_kw_spec(omg,binedge_kappa,KE_kappaT_eff,Ndays,spec_name,casename,title_str,clim_min,clim_max,save_img_path):
    omg_eff=omg[0:Nt//2]
    kappa_1d=0.5*(binedge_kappa[0:-1]+binedge_kappa[1:])*1000
    fig, ax = plt.subplots()
    plt.yscale('log'),plt.ylim([0.06,4])
    plt.xscale('log'),plt.xlim([0.004,kappa_1d.max()/np.sqrt(2.0)])
    K,W=np.meshgrid(kappa_1d,omg_eff)
    plt.pcolormesh(kappa_1d,omg_eff,tch.log(tch.tensor(W)*tch.tensor(K)*KE_kappaT_eff)/tch.log(tch.tensor(10.0)),cmap='gist_ncar')
    plt.colorbar()
    plt.clim(clim_min,clim_max)
    plt.title(title_str)
    plt.ylabel('cpd')
    plt.xlabel('cycle per km')
    fig.savefig(save_img_path+'/'+casename+'_'+loc_str+'_'+str(Ndays)+'days_'+spec_name+'_spec_kw_Hann.png')

# Uvel__rel_uvar_9__0000700000_0000720446_file_001_of_002_014days_layer_1.pt

args_cls = parser_diag_cls.parse_args()

dict_load_ts={
    'dx' : int(args_cls.dx),
    'dy' : int(args_cls.dy),
    'Lx' : int(args_cls.Lx), 
    'Ly' : int(args_cls.Ly),
    'path_scratch' : args_cls.path_scratch,
    'path_results' : args_cls.path_results,
    'groupname' : args_cls.groupname,
    'casename' : args_cls.casename,
    'start_ind' : args_cls.start_ind,
    'end_ind' : args_cls.end_ind,
    'tape_days' : args_cls.tape_days, # tape every xx days
    'Fs' : args_cls.Fs, # samples a day
    'dt_model' : args_cls.dt_model, # [sec] timestep in the simulations
    'z_target' : np.array(args_cls.z_target.split(','),dtype=int), # [1] for surface
    'name_fields' : args_cls.name_fields.split(',') # name used as the head of saved fields (tensor)
    }
print('Parameters read from at runtime')
print(dict_load_ts)

#%%
dict_load_ts={
    'dx' : 500.0,
    'dy' : 500.0,
    'Lx' : 480e3, 
    'Ly' : 960e3,
    'ind_z' : [1,2,3,4,5,20,60],
    'fhead' : 'Diag_snaps_UV',
    'path_scratch' : home_dir+'/scratch/MITgcm_cases/',
    'path_results' : home_dir+'/postproc/results/',
    'groupname' : '',
    'casename' : 'rel_uvar_9',
    'start_ind' : 700000,
    'end_ind' : 715000,
    'Fs' : 8, # samples a day
    'dt_model' : 75, # [sec] timestep in the simulations
    'z_target' : [1], # [1] for surface
    'name_fields' : ['Uvel', 'Vvel'] # name used as the head of saved fields (tensor)
    }
#%%
path_saved_obj=dict_load_ts['path_results']+'/saved_obj/'+dict_load_ts['groupname']\
    +'/'+dict_load_ts['casename']+'/'
len_fields=len(dict_load_ts['name_fields'])


#%%
UVvel_obj= load_obj(path_saved_obj, 'Uvel');

Lx=UV_obj.UVobj.Lx
Ly=UV_obj.UVobj.Ly
Fs=UV_obj.UVobj.Fs
Nt=len(UV_obj.UV_field[0])
Lt=1/8*Nt
Ny=1920 
Nx=960;

UV_obj=UV_from_ts(UVvel_obj,dict_load_ts['path_results'])
UV_obj.load_stitch_tensor([1,5])


#%%
plt_iz_target=0
if plt_iz_target==0:
    loc_str='surf'
elif plt_iz_target==1:
    loc_str='ssurf'
elif plt_iz_target==5:
    loc_str='sekm'
elif plt_iz_target==7:
    loc_str='cnt'
#%%
div2D=tch.tensor(ppf.get_divergence(UV_obj.UV_field[0][:,plt_iz_target,:,:],UV_obj.UV_field[1][:,plt_iz_target,:,:],UV_obj.UVobj.dx,UV_obj.UVobj.dy))

# FFT
omg,binedge_kappa,KE_kappaT_eff=fcn.plot_kw_spec(div2D,Nx,Ny,Nt,Lx,Ly,Lt,Fs,bins=400,clim_min=-3,clim_max=2)
# plotting 
Ndays=Nt/Fs
title_str=r'$\kappa\omega$DIV,$\Delta x$=$\Delta y$=500m'

dict_plt_kw_spectra={
    'omg':omg,
    'binedge_kappa':binedge_kappa,
    'KE_kappaT_eff':KE_kappaT_eff,
    'Ndays':Ndays,
    'spec_name':'DIV',
    'casename':dict_load_ts['casename'],
    'title_str':title_str,
    'clim_min':clim_min,
    'clim_max':clim_max,
    'save_img_path':save_img_path}

plot_kw_spec(**dict_plt_kw_spectra)

#%% Vorticity
div2D=tch.tensor(ppf.get_divergence(UV_obj.UV_field[0][:,plt_iz_target,:,:],UV_obj.UV_field[1][:,plt_iz_target,:,:],UV_obj.UVobj.dx,UV_obj.UVobj.dy))

# FFT
omg,binedge_kappa,KE_kappaT_eff=fcn.plot_kw_spec(div2D,Nx,Ny,Nt,Lx,Ly,Lt,Fs,bins=400,clim_min=-3,clim_max=2)
# plotting 
Ndays=Nt/Fs
title_str=r'$\kappa\omega$DIV,$\Delta x$=$\Delta y$=500m'

dict_plt_kw_spectra={
    'omg':omg,
    'binedge_kappa':binedge_kappa,
    'KE_kappaT_eff':KE_kappaT_eff,
    'Ndays':Ndays,
    'spec_name':'DIV',
    'casename':dict_load_ts['casename'],
    'title_str':title_str,
    'clim_min':clim_min,
    'clim_max':clim_max,
    'save_img_path':save_img_path}

plot_kw_spec(**dict_plt_kw_spectra)

