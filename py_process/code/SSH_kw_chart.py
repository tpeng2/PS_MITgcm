#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:15:38 2021
Make wave-speactra for 
@author: tpeng2
"""

#%% preset functions from postprocess_fcns.py
import sys
sys.path.append('/home/tpeng/Dropbox/MITgcm/PyPostproc/Tests')
from PPFCN import *
import PPFCN.ppfcns as ppf
# from PPFCN.ppfcns import *

#%% define functions
tch.pi=tch.tensor(np.pi)
if tch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  
#%% model parameter
# Nx=480; Ny=960; Nz=100;
Nx=960; Ny=1920; Nz=100;
Lx = 480e3; Ly = 960e3; Lz = 4000;
# Coriolis 
f_0=-1e-4;
beta=1.5e-11;
dy=Ly/Ny;
y_u=tch.linspace(dy/2,Ny*dy-dy/2,Ny)
f=f_0+beta*y_u
y_plot=15
T_iner=tch.abs(2*tch.pi/f[y_plot]/86400) # day at y=480km
cpd_iner=1/T_iner
radius_iner=1.5/f[y_plot] # [m] at y=480km
#%%
# filepath='/scratch/MITgcm_cases/woce1km_run/';
# filepath='/Users/tpeng2/scratch/MITgcm_cases/woce500_720c/AGU/';
filepath='/data/MITgcm_cases/exf_uvar/AGU/';
fname=sorted(glob.glob(filepath+'SSH*data',recursive=True));
save_img_path=home_dir+'/postproc/img/'
if not os.path.exists(save_img_path):
        os.makedirs(save_img_path)
#%%
# file index: filename[i][-15:-5]
start_ind= 315150;
end_ind= 393334;

filename=ppf.trim_fname(fname,'data',start_ind,end_ind)
print(str(len(filename)//8)+' days')
#%%

# filename=filename[:]
# start_ind=809280
# end_ind=981840

inc_ind=144;
dt=75*inc_ind;
dt_day=dt/86400;
Fs=1/dt_day;
fid=np.arange(start_ind,end_ind,inc_ind).astype(int)
Nt=len(filename)
Lt=dt_day*Nt;
t=tch.linspace(dt_day,Lt,Nt)
# time tapering filter
Tfilter=0.5*(1-tch.cos(2*np.pi/Lt*t));
# Tfilter=tch.tensor(np.blackman(Nt));
# Tfilter=tch.tensor(signal.windows.flattop(Nt));

# initialization
SSH_orig=tch.zeros([Ny,Nx,Nt])
#%% Read SSH
print('Reading SSH file from ' +filename[1][-15:-5] +' to '+filename[-1][-15:-5] )
for i in range(Nt):
    # fileid=fid[i]
    # filename=filepath+'SSH.'+format(fileid,'010d')
    # print(filename)
    if np.mod(i,25)==0:
        print(str(np.int((i+1)/Nt*100))+'%',end=' ... ')
    elif i==Nt-1:
        print('Completed') 
    SSH_orig[:,:,i]=tch.tensor(utils.rdmds(filename[i][:-5]))
# SSH_orig=SSH.clone().detach()
SSH_mrr=tch.cat((SSH_orig,tch.flip(SSH_orig,[0])),0)
#%% stiches toghether
dSSH1=(SSH_orig[-2,:,:]-SSH_orig[-3,:,:])/dy
dSSH2=(SSH_orig[-3,:,:]-SSH_orig[-4,:,:])/dy
SSH_southfill=SSH_orig[-2,:,:]+3/4*dSSH1+1/4*dSSH2
SSH_mrr[1919,:,:]=SSH_southfill
SSH_mrr[1920,:,:]=SSH_southfill
SSH_mrr[0,:,:]=SSH_mrr[1,:,:]
SSH_mrr[-1,:,:]=SSH_mrr[-2,:,:]

# Apply filters
# SSH=SSH_mrr.clone().detach()
SSH=SSH_mrr


#%% Add filter on two ends

filter_ssh=ppf.gen_filter_1d(SSH)
for i in np.arange(len(filter_ssh)):
    filter_ssh[i]=tch.tensor(filter_ssh[i])

# SSH_filtered=SSH.clone().detach()
# SSH_filtered=SSH*filter_ssh[0][:,None,None]
SSH_filtered=SSH_mrr*filter_ssh[2][None,None,:]    

#%% FFT
fnyq=Fs/2;
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


#%%
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
hist_kappa,binedge_kappa=np.histogram(kappa[~tch.isnan(kappa)],bins=280)
dbin_kappa=np.mean(np.diff(binedge_kappa))
ind_kappa=np.arange(0,len(binedge_kappa))
# mapping to the labels (ind_kappa)
label_kappa=tch.zeros(kappa.shape)
for i in np.arange(len(ind_kappa)):
    label_kappa[kappa>binedge_kappa[i]]=i
plt.imshow(label_kappa,aspect=1/4)

#%% FFT2
# A_fft_kl=tchfft.rfftn(A,dim=[1,0])
# A_fft_klw=tchfft.fftn(SSH,dim=[0,1,2])
# EA_fft_klw=get_fftpower(A_fft_klw) #torch.Size([Nx, Ny//2+1, Nt])
# label_kappa_eff=label_kappa[0:len(l),:]

A_fft_kl=tchfft.rfftn(SSH,dim=[1,0])
A_fft_klw=tchfft.fft(A_fft_kl,dim=2)



del A_fft_kl, SSH_orig
# EA_fft_kl=get_fftpower(A_fft_kl) #torch.Size([Ny//2+1, Nx, Nt])
EA_fft_klw=get_fftpower(A_fft_klw) 
print('Shape of 2D spectra of SSH: '+str(EA_fft_klw.shape))
del A_fft_klw
kappa_crop=kappa_circ[0:Ny+1,:]
# EA_fft_klw[tch.where(kappa[0:Ny//2+1,:]>kappa_limit)]=tch.tensor(float('nan'))

label_kappa_eff=label_kappa[0:len(l)//2+1,:]

#%%create 2D wavenumber-time matrix
EA_kappaT=tch.zeros([len(hist_kappa),Nt])
print(EA_kappaT.shape)
# average 
nhist_kappa=tch.zeros(len(hist_kappa))
# go through each kappa label 
for i in np.arange(int(len(ind_kappa)//np.sqrt(2)-1)):
    # find label first
    lbl_y,lbl_x=tch.where(label_kappa_eff==i)
    EA_kappaT[i,:]=tch.sum(EA_fft_klw[lbl_y,lbl_x,:],dim=[0])
EA_kappaT_eff=EA_kappaT[:,0:Nt//2]

#%% set axis
omg_eff=omg[0:Nt//2]
kappa_1d=0.5*(binedge_kappa[0:-1]+binedge_kappa[1:])*1000
fig, ax = plt.subplots()
plt.yscale('log'),plt.ylim([0.05,4])
plt.xscale('log'),plt.xlim([0.004,kappa_1d.max()/np.sqrt(2.0)])
W,K=np.meshgrid(omg_eff,kappa_1d)
im=ax.pcolormesh(omg_eff,kappa_1d,tch.transpose(tch.log(tch.tensor(W)*tch.tensor(K)*EA_kappaT_eff)/tch.log(tch.tensor(10.0)),0,1),cmap='gist_ncar')
fig.colorbar(im, orientation='vertical')
# 1D spectra integrals
EA_1d_kappa=tch.nansum(EA_kappaT_eff,dim=1)
EA_1d_omg=tch.nansum(EA_kappaT_eff,dim=0)
plt.title(r'$k\omega$SSH, dx=dy=500m')
plt.ylabel('cpd')
plt.xlabel('cycle per km')
plt.savefig(save_img_path+'/'+casename+'SSH_spec_kw_Hann.png')
plt.show()
#%%
axs1=plt.subplot(1, 2, 1)
axs1.loglog(omg_eff,EA_1d_omg),plt.xlim([0.05,4])
plt.ylim([1e9,1e15])
plt.xlabel('cpd')
plt.title(r'$T=2\pi/f_0$')

axs2=plt.subplot(1, 2, 2)
# axs2.loglog(kappa_1d,tch.tensor(kappa_1d)**4*EA_1d_kappa),plt.xlim([0.001,1])
axs2.loglog(kappa_1d,EA_1d_kappa),plt.xlim([0.001,1])
plt.xlabel('cycle per km')
plt.ylim([1e2,1e18])
plt.savefig('SSH_spec_1d_Hann.png')
plt.title(r'$\kappa=\sqrt{k^2+l^2}$')

#%% PyTorch Nanmean
def nanmean(v, *args, inplace=False, **kwargs):
    if not inplace:
        v = v.clone()
    is_nan = tch.isnan(v)
    v[is_nan] = 0
    return v.sum(*args, **kwargs) / (~is_nan).float().sum(*args, **kwargs)
#%%
# EA_kappa_1d=tch.zeros([len(hist_kappa)])
# A_fft_kl=tchfft.fftn(SSH_filtered,dim=[0,1])
# EA_fft_kl_mean=get_fftpower(A_fft_kl)
# EA_fft_kl_mean=tch.nansum(EA_fft_kl,dim=2)/Nt
# nhist_kappa=tch.zeros(len(hist_kappa))
# # go through each kappa label 
# for i in np.arange(int(len(ind_kappa)//np.sqrt(2)-1)):
#     # find label first
#     lbl_y,lbl_x=tch.where(label_kappa_eff==i)
#     EA_kappa_1d[i]=tch.sum(EA_fft_kl_mean[lbl_y,lbl_x],dim=0)
# EA_kappa_bin=EA_kappa_1d[:]