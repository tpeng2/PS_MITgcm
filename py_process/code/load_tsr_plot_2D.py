#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import matplotlib
matplotlib.use('Agg')
home_dir=os.path.expanduser("~")
package_path=home_dir+'/MITgcm_post/py_process/ppmod/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
import PPFCN.fftfcn as fcn

from scipy import interpolate

#%%

import argparse
import shutil


parser = argparse.ArgumentParser(description="Read ExfSuf")
parser.add_argument('path_tsr', metavar="path_tsr", help="Path where tensors are stored)")
parser.add_argument('fname_vort_tsr', metavar="fname_vort_tsr", help="File name (vorticity before normalization)")
parser.add_argument('fname_vort_eta_tsr', metavar="fname_vort_eta_tsr", help="File name (inferred vorticity normalize by f=f_0+beta*y")
parser.add_argument('Fs', metavar="Fs", help="Sampling rate (per day)")
parser.add_argument('dt_sec', metavar="dt_sec", help="Time step size (in seconds)")
parser.add_argument('flg_plot_zeta', metavar="flg_plot_zeta", help="Plot relative vorticity Lap eta ('Y' or 'N')")
parser.add_argument('flg_plot_lapeta', metavar="flg_plot_lapeta", help="Plot inferred vorticity Lap eta ('Y' or 'N')")
parser.add_argument('flg_plot_diff', metavar="flg_plot_diff", help="Plot difference in vorticities ('Y' or 'N')")

args = parser.parse_args()

path_tsr=args.path_tsr
fname_vort_tsr=args.fname_vort_tsr
fname_vort_eta_tsr=args.fname_vort_eta_tsr
Fs=int(args.Fs)
dt_sec=int(args.dt_sec)
flg_plot_zeta=args.flg_plot_zeta
flg_plot_lapeta=args.flg_plot_lapeta
flg_plot_diff=args.flg_plot_diff

print('Read path_tsr as '+args.path_tsr)
print('Read fname_vort_tsr as '+args.fname_vort_tsr)
print('Read fname_vort_eta_tsr as '+args.fname_vort_eta_tsr)
print('Read Fs as '+args.Fs)
print('Read dt as '+args.dt_sec)
#%%

# path_tsr=home_dir+'/postproc/tsr/'
# fname_vort_tsr='Vort_surf_rel_ucon_9_0000480118_0000549958.pt'
# fname_vort_eta_tsr='Vort_eta_surf_rel_ucon_9_0000480118_0000549958.pt'
# Fs=8
# dt_sec=75
# flg_plot_zeta='N'
# flg_plot_lapeta='Y'
# flg_plot_diff='N'
#%%
if(flg_plot_zeta=='Y' or flg_plot_lapeta==1):
    print('Loading vorticity tensor: '+fname_vort_tsr)
    vort=tch.load(path_tsr+fname_vort_tsr)
    print('Size of vorticity tensor: '+str(vort.shape))
if(flg_plot_lapeta=='Y' or flg_plot_lapeta==1):
    print('Loading Lap(eta) tensor: '+str(fname_vort_eta_tsr))
    vort_eta=tch.load(path_tsr+fname_vort_eta_tsr)
    print('Size of Lap(eta) tensor: '+str(vort_eta.shape))
#%%
start_ind_str=fname_vort_tsr[-24:-14]
end_ind_str=fname_vort_tsr[-13:-3]
start_ind_eta_str=fname_vort_tsr[-24:-14]
end_ind_eta_str=fname_vort_tsr[-13:-3]
    
start_ind=int(start_ind_str)
end_ind=int(end_ind_str)

df_ind=86400/dt_sec/Fs;

dx=500
dy=500
Lx=480e3
Ly=960e3
xc=np.linspace(dx/2,Lx-dx/2,int(Lx//dx))
yc=np.linspace(dy/2,Ly-dy/2,int(Ly//dy))

f,f_2D=ppf.gen_f_coriolis(-1e-4,1.5e-11,xc,yc)
if flg_plot_zeta=='Y':
    Nfile=len(vort)
elif flg_plot_lapeta=='Y':
    Nfile=len(vort_eta)



#%% Frequency filtering
T_wind=60480
f_iner=86400/T_wind
def save_filtered_2D_tsr(M,Fs,f_cut_low,f_cut_right,path_tsr,fname_M,label_str):
    M_filtered=fcn.fourier_freq_filter(M,Fs,f_cut_low,f_cut_right)
    tch.save(M_filtered,path_tsr+label_str+fname_M,pickle_protocol=4)
    del M_filtered
    print('Filtered field: '+label_str+'_'+fname_M)

f_cut_left=0
f_cut_mid=0.1*f_iner
f_cut_right=2*f_iner

# vort_eta
save_filtered_2D_tsr(vort_eta,Fs,f_cut_left,f_cut_mid,path_tsr,fname_vort_eta_tsr,'fltd_low_')
save_filtered_2D_tsr(vort_eta,Fs,f_cut_mid,f_cut_right,path_tsr,fname_vort_eta_tsr,'fltd_mid_')
save_filtered_2D_tsr(vort_eta,Fs,f_cut_right,Fs,path_tsr,fname_vort_eta_tsr,'fltd_high_')

# vort
save_filtered_2D_tsr(vort,Fs,f_cut_left,f_cut_mid,path_tsr,fname_vort_tsr,'fltd_low_')
save_filtered_2D_tsr(vort,Fs,f_cut_mid,f_cut_right,path_tsr,fname_vort_tsr,'fltd_mid_')
save_filtered_2D_tsr(vort,Fs,f_cut_right,Fs,path_tsr,fname_vort_tsr,'fltd_high_')

path_save_img=home_dir+'/postproc/img/'

#%%
if(flg_plot_zeta=='Y' or flg_plot_lapeta=='Y' or flg_plot_diff=='Y'):
    for i in range(Nfile):
        day='{ind:07.3f}'.format(ind=i*1/Fs)
        print('Plotting '+ str(i)+' of '+ str(Nfile) + ' files.')
        if (flg_plot_zeta=='Y'):
            fpath=save_img_path=home_dir+'/postproc/img/'+fname_vort_tsr[:25]
            if not os.path.exists(fpath):
                os.makedirs(fpath)
            ppf.plot_2d_colorbar(vort[i,15:-15]/f_2D[15:-15],xc/1000,yc[15:-15]/1000,'x [km]','y [km]',
                            r'Relative vorticity $\zeta/f$','bwr',[-1,1],plt_asp=1,
                            fname=fname_vort_tsr[:25]+str(start_ind)+'_'+str(day)+'_days',ftype='.png',fpath=home_dir+'/postproc/img/')
        if (start_ind_str == start_ind_eta_str and flg_plot_lapeta=='Y'):
            fpath=save_img_path=home_dir+'/postproc/img/'+fname_vort_eta_tsr[:25]
            if not os.path.exists(fpath):
                os.makedirs(fpath)
            ppf.plot_2d_colorbar(vort_eta[i,15:-15],xc/1000,yc[15:-15]/1000,'x [km]','y [km]',
                             r'Inferred vorticity $\nabla^2 \eta /f^2 g$','bwr',[-1,1],plt_asp=1,
                             fname=fname_vort_eta_tsr[:25]+str(start_ind)+'_'+str(day)+'_days',ftype='.png',fpath=home_dir+'/postproc/img/')
            if (flg_plot_diff=='Y' or flg_plot_diff==1):
                fpath=save_img_path=home_dir+'/postproc/img/d_'+fname_vort_eta_tsr[:25]
                if not os.path.exists(fpath):
                    os.makedirs(fpath)
                ppf.plot_2d_colorbar(vort[i,15:-15]/f_2D[15:-15]-vort_eta[i,15:-15],xc/1000,yc[15:-15]/1000,'x [km]','y [km]',
                                 r'Differences: $\zeta/f-\nabla^2 \eta /f^2 g$','bwr',[-1,1],plt_asp=1,
                                 fname='d_'+fname_vort_tsr[:25]+str(start_ind)+'_'+str(day)+'_days',ftype='.png',fpath=home_dir+'/postproc/img/')
