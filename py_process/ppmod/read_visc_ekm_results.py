#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 01:07:43 2021
Read Viscosity Test Results
@author: tpeng
"""

import sys

import os
home_dir=os.path.expanduser("~")
sys.path.append(home_dir+'/data/'+'/Dropbox/MITgcm/PyPostproc/Tests')

from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate
#%%
groupname=''
casename='rel_uvar2'
# path_UV=home_dir+'/data/MITgcm/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
path_UV=home_dir+'/data/ExtSurf/'+groupname+'/'+casename+'/'

# vertical grids
RF=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/RF')
RC=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/RC')
DRC=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/DRC')
DRF=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/DRF')

Usekm_fname = ppf.search_binary(path_UV,'Usekm*')
# Index is stored as Usekm_fname[i][-10:]

start_ind= 363622;
end_ind= 363722;
Fs=8; # 8 slides per day

Usekm_fname_trimmed = ppf.trim_fname(Usekm_fname,'',start_ind,end_ind)
Nfile_UV=len(Usekm_fname_trimmed);

Nx = 960; 
Ny = 1920; 
Nz=100;

# Usekm=np.zeros((Ny,Nx,Nfile_UV))
# Vsekm=np.zeros((Ny,Nx,Nfile_UV))
# Psekm=np.zeros((Ny,Nx,Nfile_UV))
# # Usurf=np.zeros((Ny,Nx,Nfile_UV))
# # Vsurf=np.zeros((Ny,Nx,Nfile_UV))
# # Wcnt=np.zeros((Ny,Nx,Nfile_UV))
# # Tsurf=np.zeros((Ny,Nx,Nfile_UV))

# Stats_1D = np.zeros((Nfile_UV,Nz))
file_ind=[]
# for i in range(Nfile_UV):
#     file_ind.append(Usekm_fname_trimmed[i][-10:])
#     Usekm[:,:,i]=np.fromfile(path_UV+'Usekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     Vsekm[:,:,i]=np.fromfile(path_UV+'Vsekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Usurf[:,:,i]=np.fromfile(path_UV+'Usurf'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Vsurf[:,:,i]=np.fromfile(path_UV+'Vsurf'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     Psekm[:,:,i]=np.fromfile(path_UV+'Psekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Ucnt[:,:,i]=np.fromfile(path_UV+'Ucnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Vcnt[:,:,i]=np.fromfile(path_UV+'Vcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Wcnt[:,:,i]=np.fromfile(path_UV+'Wcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Wcnt[:,:,i]=np.fromfile(path_UV+'Wcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
#     # Wxzcnt[:,:,i]=np.fromfile(path_UV+'Wxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
#     # Uxzcnt[:,:,i]=np.fromfile(path_UV+'Uxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
#     # Wssurf[:,:,i]=np.fromfile(path_UV+'Wssurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)
#     # Tsurf[:,:,i]=np.fromfile(path_UV+'Tsurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)

# # vorticity
dx=500; dy = 500;
# zeta_all=ppf.get_vorticity(Usekm,Vsekm,dx,dy)
# # zeta_all_surf=ppf.get_vorticity(Usurf,Vsurf,dx,dy)
# # dRC=np.diff(RC,axis=0)
# # dUxz_dz=np.diff(Uxzcnt,axis=0)/np.diff(RC,axis=0)
# # dWxz_dx=np.diff(Wxzcnt,axis=1)/dx
# # dU_dz_ff=np.zeros((Nz,Nx,Nfile_UV))
# # dW_dx_ff=np.zeros((Nz,Nx,Nfile_UV))
# # dU_dz_ff[1:,:]=dUxz_dz
# # dW_dx_ff[:,1:]=dWxz_dx
# # dU_dz_ff[0,:]=Uxzcnt[0,:]/(dRC[0]/2)
# # dW_dx_ff[:,0]=(Wxzcnt[:,0]-Wxzcnt[:,-1])/(dx)
# # vort_xz=dU_dz_ff-dW_dx_ff
#%% 
# UsekmFFT=tchfft.fft(tch.tensor(Usekm),dim=1); 
# UFFT2=UsekmFFT*tch.conj(UsekmFFT)
# VsekmFFT=tchfft.fft(tch.tensor(Vsekm),dim=1); 
# VFFT2=VsekmFFT*tch.conj(VsekmFFT)

# # Vort_surfFFT=tchfft.fft(tch.tensor(zeta),dim=1)
# # ZETAFFT2=Vort_surfFFT*tch.conj(Vort_surfFFT)
# Nk=Nx
# n_k=Nx*dx*2/np.linspace(1,Nx,num=Nx//2+1,dtype=np.int)
# KE1d=np.array(tch.sum((UFFT2[60:-60,0:Nx//2+1]+VFFT2[60:-60,0:Nx//2+1]).real,dim=[0,2])/UFFT2.shape[2])
# fig_1dfftx=plt.figure()
# plt.loglog(1/n_k,KE1d)
# plt.xlabel(r'k $[m^{-1}]$')
# plt.ylabel('Power')
# plt.title('Ufftx*conj(Ufftx)+Vfftx*conj(Vfftx): '+casename)
# plt.savefig('./images/KExmean'+casename+'.png')
# plt.ylim([1e-3,1e9])
# # fig_1dvortfftx=plt.figure()
# # plt.loglog(tch.sum((ZETAFFT2[:,0:Nx//2+1]).real,dim=[0,2]))
# # plt.xlabel('n_k')
# # plt.ylabel('Power')
# # plt.title('Ufftx*conj(enstrophy)+Vfftx*conj(enstrophy)')
# plt.savefig('./images/ENTxmean'+casename+'.png')
#%% plotting pcolormesh
#backwards index 
day=0
for i in np.linspace(0,Nfile_UV,int(np.ceil(Nfile_UV)),dtype=int):
    file_ind.append(Usekm_fname_trimmed[i][-10:])
    Usekm=np.fromfile(path_UV+'Usekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    Vsekm=np.fromfile(path_UV+'Vsekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Usurf[:,:,i]=np.fromfile(path_UV+'Usurf'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Vsurf[:,:,i]=np.fromfile(path_UV+'Vsurf'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    Psekm=np.fromfile(path_UV+'Psekm'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Ucnt[:,:,i]=np.fromfile(path_UV+'Ucnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Vcnt[:,:,i]=np.fromfile(path_UV+'Vcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Wcnt[:,:,i]=np.fromfile(path_UV+'Wcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Wcnt[:,:,i]=np.fromfile(path_UV+'Wcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Wxzcnt[:,:,i]=np.fromfile(path_UV+'Wxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
    # Uxzcnt[:,:,i]=np.fromfile(path_UV+'Uxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
    # Wssurf[:,:,i]=np.fromfile(path_UV+'Wssurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)
    # Tsurf[:,:,i]=np.fromfile(path_UV+'Tsurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)
    
    print(i)
    zeta=ppf.get_vorticity(Usekm,Vsekm,dx,dy)
    
    zeta_at_eta=np.zeros(zeta.shape)
    
    zeta_intp=1/4*(zeta[0:-1,0:-1]+zeta[1:,1:]+zeta[0:-1,1:]+zeta[1:,0:-1])
    
    day=day+1/8
    
    Lx=480e3;
    Ly=960e3;
    dx=500;
    dy=500;
    
    # test
    # zeta_shift,residual_max,x_eta,y_eta=ppf.interp_eta_zeta(zeta,250,250,Lx,Ly,500,1e-3)  
    #% Load SSH 
    # path_SSH=home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/AGU/'
    # SSH_fname=path_SSH+'SSH.'+file_ind[i]
    # SSH_orig=utils.rdmds(SSH_fname)
    SSH_orig=Psekm/9.81
    eta_at_zeta,res_eta_max,x_zeta,y_zeta=ppf.interp_eta_zeta(SSH_orig,-250,-250,Lx,Ly,30,2e-6)  
    # #%
    Lap_eta,Str_eta=ppf.calc_Laplacian(eta_at_zeta, dx, dy)
    
    
    #vort=Lap_eta/(f[1:-1,None]*f[1:-1,None])*9.81
    # Plot field
    # coriolis
    f,f_2D=ppf.gen_f_coriolis(-1e-4, 1.5e-11, x_zeta, y_zeta)
    Vort_eta=Lap_eta/(f_2D*f_2D)*9.81
    # ppf.plot_2d_colorbar(Vort_eta[25:-25,:],x_zeta/1e3,y_zeta[25:-25]/1e3,\
    #                      'x [km]','y [km]',r'Vorticity from L($\eta^{(at~\zeta)})/f^2*g$','bwr',(-1.5,1.5))
    # ppf.plot_2d_colorbar(zeta[25:-25,:]/f_2D[25:-25],x_zeta/1e3,y_zeta[25:-25]/1e3,\
    #                      'x [km]','y [km]',r'Vorticity from $\zeta/f$','bwr',(-1.5,1.5))
    
    opt_save_Lap_eta={"fname":'Lap_sekm_eta_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}
    opt_save_zeta_ekm={"fname":'zeta_sekm_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}
    opt_save_zeta={"fname":'zeta_surf_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}
    opt_save_dzeta={"fname":'dzeta_sekm_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}
    opt_save_ekm_eta={"fname":'dekm_zeta_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}
    opt_save_w={"fname":'w_cnt_'+casename+'_'+''+file_ind[i],"fpath":'/home/tpeng/Documents/images/',"ftype":'.png'}

    ppf.plot_2d_colorbar(Vort_eta[60:-60,:],x_zeta/1e3,y_zeta[60:-60]/1e3,\
                          'x [km]','y [km]',r'Vorticity from L($\eta^{(at~\zeta)})/f^2*g$'+'; ' +casename+'#'+file_ind[i],'bwr',(-1,1),**opt_save_Lap_eta)
    
    # ppf.plot_2d_colorbar(zeta[60:-60,:]/f_2D[60:-60,:],x_zeta/1e3,y_zeta[60:-60]/1e3,\
    #                       'x [km]','y [km]',r'Vorticity from $\zeta_{ekm}/f$'+'; ' +casename+'#'+file_ind[i]+', day: +'+'{ind_disp:06.3f}'.format(ind_disp=day),'bwr',(-1.0,1.0),**opt_save_zeta_ekm)
    # ppf.plot_2d_colorbar(zeta_all_surf[60:-60,:,i]/f_2D[60:-60,:],x_zeta/1e3,y_zeta[60:-60]/1e3,\
    #                        'x [km]','y [km]',r'Vorticity from $\zeta/f$'+'; ' +casename+'#'+file_ind[i]+', day: +'+str(day),'bwr',(-1.5,1.5),**opt_save_zeta)
    ppf.plot_2d_colorbar(zeta[60:-60,:]/f_2D[60:-60,:]-Vort_eta[60:-60,:],x_zeta/1e3,y_zeta[60:-60]/1e3,\
                          'x [km]','y [km]',r'$\zeta/f$-L($\eta^{(at~\zeta)})/f^2*g$ at -70m'+'; ' +casename+'#'+file_ind[i],'bwr',(-1,1),**opt_save_dzeta)
    
    # ppf.plot_2d_colorbar((zeta_all_surf[60:-60,:,i]-zeta[60:-60,:])/f_2D[60:-60,:],x_zeta/1e3,y_zeta[60:-60]/1e3,\
    #                      'x [km]','y [km]',r'$\Delta \zeta_{ekm}$'+'; ' +casename+'#'+file_ind[i]+', day: +'+str(day),'bwr',(-1.5,1.5),**opt_save_ekm_eta)
    # ppf.plot_2d_colorbar(Wcnt[60:-60,:,i],x_zeta/1e3,y_zeta[60:-60]/1e3,\
    #                       'x [km]','y [km]',r'$w$ [m/s]'+'; ' +casename+'#'+file_ind[i],'PuOr_r',(-0.001,0.001),plt_asp=1,**opt_save_w)
    del Usekm; del Vsekm; del Psekm;
#%% Filter

from scipy import special,fftpack
import numpy as np

def circularLowpassKernel(omega_c, N):  # omega = cutoff frequency in radians (pi is max), N = horizontal size of the kernel, also its vertical size, must be odd.
  kernel = np.fromfunction(lambda x, y: omega_c*special.j1(omega_c*np.sqrt((x - (N - 1)/2)**2 + (y - (N - 1)/2)**2))/(2*np.pi*np.sqrt((x - (N - 1)/2)**2 + (y - (N - 1)/2)**2)), [N, N])
  kernel[(N - 1)//2, (N - 1)//2] = omega_c**2/(4*np.pi)
  return kernel

kernelN = 960
dxkm=0.5 #km
filter_width=10#km
nfac_filter=int(Nx/filter_width*dxkm)
omega_c = np.pi*1/nfac_filter
kernel = circularLowpassKernel(omega_c, kernelN)
kernel_sft=fftpack.fftshift(kernel)/np.max(kernel)

plt.imshow(kernel)

# Crop field
zeta_cnt=zeta[480:-480,:]
eta_at_zeta_cnt=eta_at_zeta[480:-480,:]

FT_zeta_cnt=tchfft.fftn(tch.tensor(zeta_cnt),dim=[0,1])
FT_eta_cnt=tchfft.fftn(tch.tensor(eta_at_zeta_cnt),dim=[0,1])

FT_zeta_cnt_LP=FT_zeta_cnt*kernel_sft
FT_eta_cnt_LP=FT_eta_cnt*kernel_sft

zeta_cnt_LP=tch.real(tchfft.ifftn(FT_zeta_cnt_LP,dim=[0,1]))
eta_cnt_LP=tch.real(tchfft.ifftn(FT_eta_cnt_LP,dim=[0,1]))

zeta_LP_plot=zeta_cnt_LP/f_2D[480:-480]

ppf.plot_2d_colorbar(zeta_LP_plot,x_zeta/1e3,y_zeta[480:-480]/1e3,\
                         'x [km]','y [km]',r'$\zeta_{LP,'+str(filter_width)+'km}$'+'; ' +'ref''','bwr',(-1,1),**opt_save_eta)
   
iylim=480+filter_width
Lap_eta_LP,Str_eta_LP=ppf.calc_Laplacian(eta_cnt_LP, dx, dy)
Vort_eta_LP=Lap_eta_LP/(f_2D[480:-480]*f_2D[480:-480])*9.81    

ppf.plot_2d_colorbar(tch.tensor(Vort_eta_LP[filter_width:-filter_width]),x_zeta/1e3,y_zeta[iylim:-iylim]/1e3,\
                         'x [km]','y [km]',r'$\nabla^2\eta_{LP,'+str(filter_width)+'km}$'+'; ' +'ref','bwr',(-1,1),**opt_save_eta)
dzeta_LP_plt=tch.tensor(Vort_eta_LP)-zeta_LP_plot
ppf.plot_2d_colorbar(dzeta_LP_plt[filter_width:-filter_width],x_zeta/1e3,y_zeta[iylim:-iylim]/1e3,\
                         'x [km]','y [km]',r'$\nabla^2\eta_{LP,'+str(filter_width)+'km}-\zeta_{LP,'+str(filter_width)+'km}$'+'; ' +'ref','bwr',(-1,1),**opt_save_eta)
