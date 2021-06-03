#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 13:45:46 2021

2D spatial FFT 
Spatial and temporal spectra of U and V

@author: tpeng
"""
import sys
import os
home_dir=os.path.expanduser("~")
package_path='/home/tpeng/MITgcm_post/py_process/Test_exf/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate
#%%

groupname=''
casename='rel_uvar_9'
loc_str='sekm'
# path_UV=home_dir+'/data/MITgcm/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
path_UV=home_dir+'/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
save_img_path=home_dir+'/postproc/img/'
if not os.path.exists(save_img_path):
        os.makedirs(save_img_path)

# # vertical grids
# RF=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/RF')
# RC=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/RC')
# DRC=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/DRC')
# DRF=utils.rdmds(home_dir+'/data/MITgcm_cases/'+groupname+'/'+casename+'/DRF')

Usurf_fname = ppf.search_binary(path_UV,'Usurf*')
# Index is stored as Usurf_fname[i][-10:]

start_ind= 520006;
end_ind= 574406;
Fs=8; # 8 slides per day

Usurf_fname_trimmed = ppf.trim_fname(Usurf_fname,'',start_ind,end_ind)
Nfile_UV=len(Usurf_fname_trimmed);
Nday=Nfile_UV//Fs
print(str(Nday)+' days of sample will be processed.')

Nx = 960; 
Ny = 1920; 
Nz=100;

Usurf=np.zeros((Ny,Nx,Nfile_UV))
Vsurf=np.zeros((Ny,Nx,Nfile_UV))
# Ucnt=np.zeros((Ny,Nx,Nfile_UV))
# Vcnt=np.zeros((Ny,Nx,Nfile_UV))
# Wcnt=np.zeros((Ny,Nx,Nfile_UV))
# Wxzcnt=np.zeros((Nz,Nx,Nfile_UV))
# Uxzcnt=np.zeros((Nz,Nx,Nfile_UV))
# # Wssurf=np.zeros((Ny,Nx,Nfile_UV))
# Tsurf=np.zeros((Ny,Nx,Nfile_UV))

Stats_1D = np.zeros((Nfile_UV,Nz))
file_ind=[]
for i in range(Nfile_UV):
    file_ind.append(Usurf_fname_trimmed[i][-10:])
    Usurf[:,:,i]=np.fromfile(path_UV+'U'+loc_str+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    Vsurf[:,:,i]=np.fromfile(path_UV+'V'+loc_str+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Ucnt[:,:,i]=np.fromfile(path_UV+'Ucnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Vcnt[:,:,i]=np.fromfile(path_UV+'Vcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Wcnt[:,:,i]=np.fromfile(path_UV+'Wcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    # Wxzcnt[:,:,i]=np.fromfile(path_UV+'Wxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
    # Uxzcnt[:,:,i]=np.fromfile(path_UV+'Uxzcnt'+'.'+file_ind[i],dtype=np.float32).reshape(Nz,Nx)
    # Wssurf[:,:,i]=np.fromfile(path_UV+'Wssurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)
    # Tsurf[:,:,i]=np.fromfile(path_UV+'Tsurf'+'.'+file_ind,dtype=np.float32).reshape(Ny,Nx)
print('files are loaded')
# vorticity
dx=500; dy = 500;
zeta_all=ppf.get_vorticity(Usurf,Vsurf,dx,dy)


dt=1/Fs
Usurf=tch.tensor(Usurf).detach().clone()
Vsurf=tch.tensor(Vsurf).detach().clone()

def mirror_field_in_y(M,dy):
    Nx=M.shape[1]
    Ny=M.shape[0]
    Nt=M.shape[2]
    dM1=(M[-2,:,:]-M[-3,:,:])/dy
    dM2=(M[-3,:,:]-M[-4,:,:])/dy
    M_southfill=M[-2,:,:]+3/4*dM1+1/4*dM2
    M_mrr=tch.zeros((Ny*2,Nx,Nt))
    M_mrr[:Ny,:,:]=M.clone().detach();
    M_mrr[Ny:,:,:]=tch.flipud(M);
    M_mrr[Ny-1,:,:]=M_southfill
    M_mrr[Ny,:,:]=M_southfill
    M_mrr[0,:,:]=M_mrr[1,:,:]
    M_mrr[-1,:,:]=M_mrr[-2,:,:]
    return M_mrr

dy=500 # [m]
Ny=1920 
Nx=960;



Usurf_mrr=mirror_field_in_y(Usurf,dy)
Vsurf_mrr=mirror_field_in_y(Vsurf,dy)

#%% Filter in time
U=Usurf_mrr
V=Vsurf_mrr
filter_U=ppf.gen_filter_1d(U)
filter_V=ppf.gen_filter_1d(V)
for i in np.arange(len(filter_U)):
    filter_U[i]=tch.tensor(filter_U[i])
    filter_V[i]=tch.tensor(filter_V[i])

# SSH_filtered=SSH.clone().detach()
# SSH_filtered=SSH*filter_ssh[0][:,None,None]
U_filtered=U*filter_U[2][None,None,:]    
V_filtered=V*filter_V[2][None,None,:]    
print('Filters are added to U and V fields')
#%% FFT function copied
fnyq=Fs/2;
Nt=Nfile_UV
f_slow=Fs/Nt;
# === def fun():  wavenumber array from FFT length ===
def gen_n_vec(Nt):
    if np.mod(Nt//2,2)==0: # even
        n=np.concatenate([np.linspace(0,Nt//2-1,Nt//2),np.linspace(-(Nt//2),-1,Nt//2)])
    else: #odd
        n=np.concatenate([np.linspace(0,Nt//2,Nt//2+1),np.linspace(-(Nt//2),-1,Nt//2)])
    return tch.tensor(n,dtype=int).clone().detach()
# === def fun():  temporal frequency array, scaled with the length ===
def n_to_freq(n,Fs):
    Nsmp=len(n);
    freq=n/Nsmp*Fs
    return freq.clone().detach()
# === def fun():  Get power spectrum ===
def get_fftpower(A):
    A_power=tch.abs(A)*tch.abs(A.conj())
    return A_power.clone().detach()

# === def fun(): get 2-D wavenumber map with a different aspect ratio
def kl_to_kappa(k,l):
    Nk=len(k)
    Nl=len(l)
    [kk,ll]=tch.meshgrid([k,l])
    kappa=tch.sqrt(kk**2+ll**2)
    return kappa.clone().detach()

#%%#%%
Lx=480e3
Ly=960e3
Lt=1/8*Nt
nk=gen_n_vec(Nx)
nl=gen_n_vec(2*Ny)
nomg=gen_n_vec(Nt)
# scale with its length
k=n_to_freq(nk,Nx/Lx) 
l=n_to_freq(nl,Ny/Ly)
omg=n_to_freq(nomg,Nt/Lt)
# get 2D kappa map
kappa=kl_to_kappa(l,k)
# remove corner
print('Length of k: '+str(kappa.shape[1]))
print('Length of l: '+str(kappa.shape[0]))
kappa_limit=tch.max(kappa)/np.sqrt(2)
kappa_circ=kappa.clone().detach()
kappa_circ[tch.where(kappa>kappa_limit)]=tch.tensor(float('nan'))

#%%  bin 2D kappa map
# create an index map
hist_kappa,binedge_kappa=np.histogram(kappa[~tch.isnan(kappa)],bins=400)
dbin_kappa=np.mean(np.diff(binedge_kappa))
ind_kappa=np.arange(0,len(binedge_kappa))
# mapping to the labels (ind_kappa)
label_kappa=tch.zeros(kappa.shape)
for i in np.arange(len(ind_kappa)):
    label_kappa[kappa>binedge_kappa[i]]=i
#plt.imshow(label_kappa,aspect=1/4)
print('2D kappa binned')
#%% FFT2
# A_fft_kl=tchfft.rfftn(A,dim=[1,0])
# A_fft_klw=tchfft.fftn(SSH,dim=[0,1,2])
# EA_fft_klw=get_fftpower(A_fft_klw) #torch.Size([Nx, Ny//2+1, Nt])
# label_kappa_eff=label_kappa[0:len(l),:]

U_fft_kl=tchfft.rfftn(U,dim=[1,0])
U_fft_klw=tchfft.fft(U_fft_kl,dim=2)

del U_fft_kl, U

V_fft_kl=tchfft.rfftn(V,dim=[1,0])
V_fft_klw=tchfft.fft(V_fft_kl,dim=2)

# EA_fft_kl=get_fftpower(A_fft_kl) #torch.Size([Ny//2+1, Nx, Nt])
E_U_fft_klw=get_fftpower(U_fft_klw) 
del U_fft_klw
E_V_fft_klw=get_fftpower(V_fft_klw) 
del V_fft_klw

print('Shape of 2D spectra: '+str(E_U_fft_klw.shape))
kappa_crop=kappa_circ[0:Ny+1,:]
# EA_fft_klw[tch.where(kappa[0:Ny//2+1,:]>kappa_limit)]=tch.tensor(float('nan'))

label_kappa_eff=label_kappa[0:len(l)//2+1,:]

#%%create 2D wavenumber-time matrix
KE_kappaT=tch.zeros([len(hist_kappa),Nt])
print(KE_kappaT.shape)
# average 
nhist_kappa=tch.zeros(len(hist_kappa))
# go through each kappa label 
for i in np.arange(int(len(ind_kappa)//np.sqrt(2)-1)):
    # find label first
    lbl_y,lbl_x=tch.where(label_kappa_eff==i)
    KE_kappaT[i,:]=tch.sum(E_U_fft_klw[lbl_y,lbl_x,:]+E_V_fft_klw[lbl_y,lbl_x,:],dim=[0])
KE_kappaT_eff=KE_kappaT[:,0:Nt//2]
#%% set axis
omg_eff=omg[0:Nt//2]
kappa_1d=0.5*(binedge_kappa[0:-1]+binedge_kappa[1:])*1000
fig, ax = plt.subplots()
plt.yscale('log'),plt.ylim([0.06,4])
plt.xscale('log'),plt.xlim([0.004,kappa_1d.max()/np.sqrt(2.0)])
W,K=np.meshgrid(omg_eff,kappa_1d)
plt.pcolormesh(kappa_1d,omg_eff,tch.transpose(tch.log(tch.tensor(W)*tch.tensor(K)*KE_kappaT_eff)/tch.log(tch.tensor(10.0)),0,1),cmap='gist_ncar')
plt.colorbar()
plt.clim(2,14)

plt.title(r'$\kappa\omega$KE,$\Delta x$=$\Delta y$=500m')
plt.ylabel('cpd')
plt.xlabel('cycle per km')
plt.savefig(save_img_path+'/'+casename+'_'+loc_str+'_'+str(Nday)+'days_KE_spec_kw_Hann.png')
# plt.show()

#%%
# 1D spectra integrals
KE_1d_kappa=tch.nansum(tch.tensor(K)*KE_kappaT_eff,dim=1)/tch.nansum(tch.tensor(kappa_1d))
EA_1d_omg=tch.nansum(tch.tensor(W)*KE_kappaT_eff,dim=0)/tch.nansum(tch.tensor(omg_eff))
plt.figure()
axs1=plt.subplot(1, 2, 1)
axs1.loglog(omg_eff,EA_1d_omg),plt.xlim([0.05,4])
plt.ylim([1e8,1e15])
plt.xlabel('cpd')
plt.ylabel(r'$ S$')
plt.title(r'$T=2\pi/f_0$')

axs2=plt.subplot(1, 2, 2)
axs2.loglog(kappa_1d,KE_1d_kappa),plt.xlim([0.002,1])
plt.xlabel('cycle per km')
plt.ylim([1e6,1e13])
plt.savefig('SSH_spec_1d_Hann.png')
plt.ylabel(r'$ S$')
plt.title(r'$\kappa=\sqrt{k^2+l^2}$')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(save_img_path+'/'+casename+'_'+str(Nday)+'days_KE_1Dspec_kw_Hann.png')
