#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:57:10 2021

read full-delpth 3D fields

* search full fields files
* get mean stratification
* get mean KE(z)
    -- full field
    -- +-25% y of the central channel
* calculate stratification and frequency


@author: tpeng
"""


import sys
import os
home_dir=os.path.expanduser("~")
package_path=home_dir+'/MITgcm_post/py_process/Test_exf/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate
#%%
import argparse
import shutil


parser = argparse.ArgumentParser(description="Read ExfSuf")
parser.add_argument('groupname', metavar="GROUPNAME", \
                    help="Group name (groupname)")
parser.add_argument('casename', metavar="CASENAME", \
                    help="Case name (casename)")
parser.add_argument('start_ind', metavar="START_IND", \
                    help="Start index (start_ind)")
parser.add_argument('end_ind', metavar="END_IND", \
                    help="Start index (end_ind)")
parser.add_argument('loc_str', metavar="LOC_STR", \
                    help="Location (e.g.,  'surf', 'sekm', or 'cnt'.")


args = parser.parse_args()


groupname=args.groupname
casename=args.casename
start_ind=int(args.start_ind)
end_ind=int(args.end_ind)
loc_str=args.loc_str

#%%
# groupname=''
# casename='rel_uvar_9'
# 
# start_ind= 470150;
# end_ind= 533334;
# loc_str='surf'
#
#%%
print("Case name "+casename+" is read.")
ext_savepath=home_dir+'/postproc/results/ExtSurf/'

# path_UV=home_dir+'/data/MITgcm/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
path_UV=home_dir+'/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
print("Full path of case extraction: "+path_UV)

path_rawfiles=home_dir+'/scratch/MITgcm_cases/'+groupname+'/'+casename+'/'
print("Raw file is stored at: " + path_rawfiles)

path_Diagfiles=home_dir+'/scratch/MITgcm_cases/'+groupname+'/'+casename+'/DNG/'
print("Diagnostics file is stored at: " + path_Diagfiles)



save_img_path=home_dir+'/postproc/img/'
save_tsr_path=home_dir+'/postproc/tsr/'
if not os.path.exists(save_img_path):
        os.makedirs(save_img_path)
if not os.path.exists(save_tsr_path):
        os.makedirs(save_tsr_path)
#%%
# XC,YC,XG,YG=ppf.load_horizontal_grids(path_rawfiles);
RC=utils.rdmds(path_rawfiles+'RC')
hFacC=utils.rdmds(path_rawfiles+'/hFacC');
# Load half cells
hFacC_nan=np.copy(hFacC)
hFacC_nan[np.where(hFacC==0)]=float('nan')
#%% search and trim files
T_fname_raw=ppf.search_mdsdata(path_rawfiles,'T');
T_fname_trimmed=ppf.trim_fname(T_fname_raw,'.data',start_ind,end_ind)

Diag_UV_fname=ppf.search_mdsdata(path_Diagfiles,'Diag*UV*');
Diag_UV_trimmed=ppf.trim_fname(T_fname_raw,'.data',start_ind,end_ind)

#%% Loop on raw files
Nfile_raw=len(T_fname_trimmed)
mean_T=np.zeros([len(RC),Nfile_raw])
rawfile_ind=np.zeros(Nfile_raw)
rawfind_str=[]
for i in range(Nfile_raw):
    rawfile_ind[i]=int(T_fname_trimmed[i][-10-len('.data'):-len('.data')])
    rawfind_str.extend(["{f_ind:010d}".format(f_ind=int(rawfile_ind[i]))])
    T=utils.rdmds(T_fname_trimmed[i][:-len('.data')])
    mean_T[:,i]=np.nanmean(T*hFacC_nan,axis=(1,2))
    plt.plot(mean_T[:,i],RC[:,0,0])
    del T;
plt.plot(np.mean(mean_T,axis=1),RC[:,0,0],ppf.T_ref_initial,RC[:,0,0])
plt.xlabel('Temperature [deg C]'); 
plt.ylabel('z [m]')
rawfind_str.extend(['Mean'])
plt.plot(ppf.T_ref_initial,RC[:,0,0]);
rawfind_str.extend(['Initial'])
plt.legend(rawfind_str)
# save figure
start_ind_str="{f_ind:010d}".format(f_ind=int(rawfile_ind[0]))
end_ind_str="{f_ind:010d}".format(f_ind=int(rawfile_ind[-1]))
plt.savefig(save_img_path+'T1Dz._'+casename+'_'\
            +start_ind_str+'_'+end_ind_str+'.png',dpi=300)
#% Save mean_T here
tch.save(tch.tensor(mean_T),save_tsr_path+'mean_T_'+
         start_ind_str+'_'+end_ind_str+'.pt',pickle_protocol=4)
#%% Loop on Diagnostics
Nfile_Diag=len(Diag_UV_trimmed)
Nfile_tape=8
mod_tape=mod(Nfile_Diag,Nfile_tape)
if mod_tape==0:
    nfile=Nfile_Diag//Nfile_tape
    slides_per_output=np.ones(nfile)*Nfile_tape
else:
    nfile=Nfile_Diag//Nfile_tape+1
    slides_per_output=np.ones(nfile)*Nfile_tape
    slides_per_output[-1]=mod_tape

    
#%% Process U, V and KE here
