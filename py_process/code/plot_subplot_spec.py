#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 12:29:37 2021
Plotting from .pkl
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
    }
print('Parameters read from at runtime')
print(dict_load_ts)

#%%
#%%
Ndays=dict_load_ts
loc_str=['surf','sekm','-200m']
loc_str_w=['ssurf','sekm','-200m']
class spec_3locs():
    def __init__(self,casename,path_spec,loc_str):
        self.casename=casename
        self.path_spec=path_spec
        self.loc_str=loc_str
        self.nsubplot=len(loc_str)
    def find_spectrum(self,specname):
        self.specname=specname
        self.dict_load_ts=['']*self.nsubplot
        self.path_spec_file=['']*self.nsubplot 
        for i in arange(self.nsubplot):
            tmp_path_spec_file=path_spec+'/'+specname+'_'+self.loc_str[i]
            self.path_spec_file[i]=glob.glob(tmp_path_spec_file+'*.pkl')[-1]
            fileobj=open(self.path_spec_file[i],'rb')
            self.dict_load_ts[i]=pickle.load(fileobj)
        self.Ndays=self.dict_load_ts[-1]['Ndays']
    def plot_subplots(self,opt_save=1):
        fig = plt.figure(figsize=(12,8));
        fig.suptitle(r'$\kappa-\omega$ spectra of '+self.specname+'\n'+' Case: '+self.casename+', '+str(int(self.Ndays))+' days',y=0.8)
        for i in arange(self.nsubplot):
            save_img_path=self.dict_load_ts[i]['save_img_path']
            binedge_kappa=self.dict_load_ts[i]['binedge_kappa']
            omg_eff=self.dict_load_ts[i]['omg_eff']
            KE_kappaT_eff=self.dict_load_ts[i]['KE_kappaT_eff']
            clim_min=self.dict_load_ts[i]['clim_min']
            clim_max=self.dict_load_ts[i]['clim_max']
            k_power=self.dict_load_ts[i]['k_power']
            Ndays=self.dict_load_ts[i]['Ndays']
            kappa_1d=0.5*(binedge_kappa[0:-1]+binedge_kappa[1:])*1000
            K,W=np.meshgrid(kappa_1d,omg_eff)
            fig.add_subplot(1,self.nsubplot,i+1);
            plt.yscale('log'),plt.ylim([0.06,4])
            plt.xscale('log'),plt.xlim([0.004,kappa_1d.max()/np.sqrt(2.0)])
            plt.pcolormesh(K,W,tch.log(tch.tensor(W)*tch.tensor(K**k_power)*KE_kappaT_eff)/tch.log(tch.tensor(10.0)),cmap='gist_ncar')
            plt.colorbar(orientation="horizontal")
            plt.clim(clim_min,clim_max)
            plt.title(self.loc_str[i])
            plt.gca().set_aspect('equal')
            plt.ylabel(r'$\omega$')
            plt.xlabel(r'$\kappa$')
            plt.tight_layout()
        plt.savefig(save_img_path+'/'+self.specname+'_'+self.casename+'_'+str(int(Ndays))+'days_3locs'+'_2Dspec_kw_Hann.png',dpi=600)
#%%
#%%
casename=dict_load_ts['casename']
path_spec='/aos/home/tpeng/postproc/img/'
def plot_subplot(casename,specname,path_spec,loc_str):
    KE=spec_3locs(casename, path_spec, loc_str)
    KE.find_spectrum(specname)
    KE.plot_subplots()

plot_subplot(casename,'KE',path_spec,loc_str)
plot_subplot(casename,'DIV',path_spec,loc_str)
plot_subplot(casename,'ZETA',path_spec,loc_str)
plot_subplot(casename,'W',path_spec,loc_str_w)
plot_subplot(casename,'Lap_P',path_spec,loc_str)

