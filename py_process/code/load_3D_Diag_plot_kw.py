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

import argparsedict_load_ts['casename']
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
    'z_target' : np.array(args_cls.z_target.split(','),dtype=int), # [1] for surface
    'name_fields' : args_cls.name_fields.split(',') # name used as the head of saved fields (tensor)
    }
print('Parameters read from at runtime')
print(dict_load_ts)

#%%

def plot_kw_spec(omg_eff,binedge_kappa,KE_kappaT_eff,Ndays,spec_name,casename,title_str,clim_min,clim_max,k_power,loc_str,save_img_path):
    kappa_1d=0.5*(binedge_kappa[0:-1]+binedge_kappa[1:])*1000
    fig, ax = plt.subplots()
    plt.yscale('log'),plt.ylim([0.06,4])
    plt.xscale('log'),plt.xlim([0.004,kappa_1d.max()/np.sqrt(2.0)])
    K,W=np.meshgrid(kappa_1d,omg_eff)
    plt.pcolormesh(K,W,tch.log(tch.tensor(W)*tch.tensor(K**k_power)*KE_kappaT_eff)/tch.log(tch.tensor(10.0)),cmap='gist_ncar')
    plt.colorbar()
    plt.clim(clim_min,clim_max)
    plt.title(title_str)
    plt.ylabel(r'$\omega$')
    plt.xlabel(r'$\kappa$')
    path_save_plot=save_img_path+'/'+casename+'_'+loc_str+'_'+str(Ndays)+'days_'+spec_name+'_spec_kw_Hann.png'
    plt.savefig(path_save_plot)
    print('Plots are saved to: '+path_save_plot)
    # 1D spectra integrals
    KE_1d_kappa=tch.nansum(tch.tensor(K)*KE_kappaT_eff,dim=0)/tch.nansum(tch.tensor(kappa_1d))
    EA_1d_omg=tch.nansum(tch.tensor(W)*KE_kappaT_eff,dim=1)/tch.nansum(tch.tensor(omg_eff))
    plt.figure()
    axs1=plt.subplot(1, 2, 1)
    axs1.loglog(omg_eff,EA_1d_omg),plt.xlim([0.05,4])
    # plt.ylim([10**clim_min/np.min(K),10**clim_max/np.max(K)])
    plt.xlabel('cpd')
    plt.ylabel(r'$ S$')
    title_str_1d='_'.join((spec_name,loc_str))
    plt.title(r'$T=2\pi/f_0$, '+title_str_1d)
    axs2=plt.subplot(1, 2, 2)
    axs2.loglog(kappa_1d,KE_1d_kappa),plt.xlim([0.002,1])
    plt.xlabel('cycle per km')
    # plt.ylim([10**clim_min/np.min(W[1:]),10**clim_max/np.max(W[1:])])
    plt.savefig('SSH_spec_1d_Hann.png') 
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$ S$')
    plt.title(r'$\kappa=\sqrt{k^2+l^2}$, '+title_str_1d)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(save_img_path+'/'+casename+'_'+loc_str+'_'+str(Ndays)+'days_KE_1Dspec_kw_Hann.png')
#%%
def get_kw_plot(M_obj,M_field,clim_min,clim_max,spec_name,opt_mirror,k_power,loc_str,save_img_path,opt_plot):
    Nx=M_obj.UV_obj.nx
    Ny=M_obj.UV_obj.ny
    Nt=M_obj.nfile_total
    Lx=M_obj.UV_obj.Lx
    Ly=M_obj.UV_obj.Ly
    Fs=M_obj.UV_obj.Fs
    Lt=Nt/Fs
    omg,binedge_kappa,KE_kappaT_eff=fcn.proc_kw_spec(M_field,Nx,Ny,Nt,Lx,Ly,Lt,Fs,bins=400,opt_mirror=opt_mirror)
    omg_eff=omg[0:Nt//2]
    Ndays=Nt/Fs
    title_str=r'$\kappa\omega$'+spec_name+', '+loc_str+',$\Delta x$=$\Delta y$=500m'
    dict_plt_kw_spectra={
        'omg_eff':omg_eff,
        'binedge_kappa':binedge_kappa,
        'KE_kappaT_eff':KE_kappaT_eff,
        'Ndays':Ndays,
        'spec_name':spec_name,
        'casename':dict_load_ts['casename'],
        'title_str':title_str,
        'clim_min':clim_min,
        'clim_max':clim_max,
        'k_power':k_power,
        'loc_str':loc_str,
        'save_img_path':save_img_path}
    #% save KE matrix
    fileobj = open(save_img_path+'_'.join([spec_name,loc_str,dict_plt_kw_spectra['casename'],str(Ndays)])+'.pkl', 'wb')
    pickle.dump(dict_plt_kw_spectra, fileobj)
    fileobj.close()
    #% if plot
    if opt_plot==1:
        plot_kw_spec(**dict_plt_kw_spectra)
    return omg,binedge_kappa,KE_kappaT_eff


#%%
# dict_load_ts={
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
#     'end_ind' : 734421,
#     'Fs' : 8, # samples a day
#     'dt_model' : 75, # [sec] timestep in the simulations
#     'z_target' : [1], # [1] for surface
#     'name_fields' : ['Uvel', 'Vvel'] # name used as the head of saved fields (tensor)
#     }
# # 
path_saved_obj=dict_load_ts['path_results']+'/saved_obj/'+dict_load_ts['groupname']\
    +'/'+dict_load_ts['casename']+'/'
len_fields=len(dict_load_ts['name_fields'])


#%%
# load class structures from files

UV_obj=UV_from_ts(path_saved_obj,'Uvel')
UV_obj.load_stitch_tensor([0,4,5])

#%% KE
def plot_KE_from_obj(UV_obj,plt_iz_target,loc_str,save_img_path):
    u2d=UV_obj.field2D[0][:,plt_iz_target,:,:]
    spec_name='U'; 
    omg,binedge_kappa,KE_kappaT_eff_U=get_kw_plot(UV_obj,u2d,clim_min=2,clim_max=12,spec_name=spec_name,opt_mirror=-1,k_power=1,loc_str=loc_str,save_img_path=save_img_path,opt_plot=1)
    del u2d
    spec_name='V'; 
    v2d=u2d=UV_obj.field2D[1][:,plt_iz_target,:,:]
    omg,binedge_kappa,KE_kappaT_eff_V=get_kw_plot(UV_obj,u2d,clim_min=2,clim_max=12,spec_name=spec_name,opt_mirror=-1,k_power=1,loc_str=loc_str,save_img_path=save_img_path,opt_plot=1)
    del v2d
    Nt=UV_obj.nfile_total
    Fs=UV_obj.UV_obj.Fs
    Lt=Nt/Fs
    omg_eff=omg[0:Nt//2]
    Ndays=Nt/Fs
    spec_name='KE'; 
    title_str=r'$\kappa\omega$'+spec_name+', '+loc_str+',$\Delta x$=$\Delta y$=500m'
    dict_plt_kw_spectra={
        'omg_eff':omg_eff,
        'binedge_kappa':binedge_kappa,
        'KE_kappaT_eff':1/2*(KE_kappaT_eff_U+KE_kappaT_eff_V),
        'Ndays':Ndays,
        'spec_name':spec_name,
        'casename':dict_load_ts['casename'],
        'title_str':title_str,
        'clim_min':4,
        'clim_max':12,
        'k_power':1,
        'loc_str':loc_str,
        'save_img_path':save_img_path}
    #% save KE matrix
    fileobj = open(save_img_path+'_'.join([spec_name,loc_str,dict_plt_kw_spectra['casename'],str(int(Ndays)),'days'])+'.pkl', 'wb')
    pickle.dump(dict_plt_kw_spectra, fileobj)
    fileobj.close()
    plot_kw_spec(**dict_plt_kw_spectra)
    
save_img_path=home_dir+'/postproc/img/'

plt_iz_target=0; loc_str='surf'
plot_KE_from_obj(UV_obj,plt_iz_target,loc_str,save_img_path)

plt_iz_target=1; loc_str='sekm'
plot_KE_from_obj(UV_obj,plt_iz_target,loc_str,save_img_path)

plt_iz_target=2; loc_str='-200m'
plot_KE_from_obj(UV_obj,plt_iz_target,loc_str,save_img_path)
#%% Divergence
plt_iz_target=0; loc_str='surf'
div2D=tch.tensor(ppf.get_divergence(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,div2D,clim_min=-3,clim_max=4,spec_name='DIV',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=1)
del div2D
#%
plt_iz_target=1; loc_str='sekm'
div2D=tch.tensor(ppf.get_divergence(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,div2D,clim_min=-3,clim_max=4,spec_name='DIV',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
del div2D

plt_iz_target=2; loc_str='-200m'
div2D=tch.tensor(ppf.get_divergence(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,div2D,clim_min=-3,clim_max=4,spec_name='DIV',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
del div2D

#%% Vorticity
plt_iz_target=0; loc_str='surf'
zeta2D=tch.tensor(ppf.get_vorticity(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,zeta2D,clim_min=-3,clim_max=4,spec_name='ZETA',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
#%
plt_iz_target=1; loc_str='sekm'
zeta2D=tch.tensor(ppf.get_vorticity(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,zeta2D,clim_min=-3,clim_max=4,spec_name='ZETA',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
plt_iz_target=2; loc_str='-200m'
zeta2D=tch.tensor(ppf.get_vorticity(UV_obj.field2D[0][:,plt_iz_target,:,:],UV_obj.field2D[1][:,plt_iz_target,:,:],UV_obj.UV_obj.dx,UV_obj.UV_obj.dy))
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(UV_obj,zeta2D,clim_min=-3,clim_max=4,spec_name='ZETA',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
del zeta2D


#%%============
# load other files
def plot_P_zeta_geo_ageo(P_obj,plt_iz_target,loc_str):
    P=P_obj.field2D[0][:,plt_iz_target,:,:]/9.81
    omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(P_obj,tch.tensor(P),clim_min=-5,clim_max=5,spec_name='Lap_P',opt_mirror=0,k_power=5,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
    # Geostrophic velocity
    ugeo,vgeo=ppf.calc_geostrophic_velocity(P,dx=500,dy=500,f_2D=f_2D)
    zeta2D_geo=tch.tensor(ppf.get_vorticity(ugeo,vgeo,P_obj.UV_obj.dx,P_obj.UV_obj.dy))
    # Ageostrophic velocity
    uageo=UV_obj.field2D[0][:,plt_iz_target,:,:]-ugeo
    vageo=UV_obj.field2D[1][:,plt_iz_target,:,:]-vgeo
    del ugeo; del vgeo
    omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(P_obj,zeta2D_geo,clim_min=-3,clim_max=4,spec_name='ZETA_geo',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
    del zeta2D_geo; del P
    zeta2D_ageo=tch.tensor(ppf.get_vorticity(uageo,vageo,P_obj.UV_obj.dx,P_obj.UV_obj.dy))
    omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(P_obj,zeta2D_ageo,clim_min=-3,clim_max=4,spec_name='ZETA_ageo',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
    del uageo; del vageo


P_obj=UV_from_ts(path_saved_obj,'P')
P_obj.load_stitch_tensor([0,4,5])
f,f_2D=ppf.gen_f_coriolis(-1.0e-4,1.5e-11,P_obj.UV_obj.xc,P_obj.UV_obj.yc)



plt_iz_target=0; loc_str='surf'
plot_P_zeta_geo_ageo(P_obj,plt_iz_target,loc_str)


plt_iz_target=1; loc_str='sekm'
plot_P_zeta_geo_ageo(P_obj,plt_iz_target,loc_str)


plt_iz_target=2; loc_str='-200m'
plot_P_zeta_geo_ageo(P_obj,plt_iz_target,loc_str)

del P_obj
del UV_obj
#%% load other files
W_obj=UV_from_ts(path_saved_obj,'W')
W_obj.load_stitch_tensor([1,4,5])


plt_iz_target=0; loc_str='ssurf'
W2D=tch.tensor(W_obj.field2D[0][:,plt_iz_target,:,:])
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(W_obj,W2D,clim_min=-3,clim_max=8,spec_name='W',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)

plt_iz_target=1; loc_str='sekm'
W2D=tch.tensor(W_obj.field2D[0][:,plt_iz_target,:,:])
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(W_obj,W2D,clim_min=-3,clim_max=8,spec_name='W',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)

plt_iz_target=2; loc_str='-200m'
W2D=tch.tensor(W_obj.field2D[0][:,plt_iz_target,:,:])
omg,binedge_kappa,KE_kappaT_eff=get_kw_plot(W_obj,W2D,clim_min=-3,clim_max=8,spec_name='W',opt_mirror=0,k_power=1,loc_str=loc_str,save_img_path=home_dir+'/postproc/img/',opt_plot=0)
del W_obj
