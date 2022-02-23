#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 15:07:04 2022

@author: tpeng
"""

# read procesed pkl file and compare cross cases
#%% Loading modules
from json import load
import os
import sys
home_dir=os.path.expanduser("~")
package_path=home_dir+'/MITgcm_post/py_process/ppmod/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate

#%%
import pickle
def load_pickle(fpath,handle):
  fobj=open(fpath,handle)
  fdata=pickle.load(fobj)
  return fdata


class spectraPKL():
  def __init__(self,casename,tDuration,pklPath):
    self.casename = casename
    self.tDuration = tDuration
    if pklPath[-1] != '/':
        pklPath += '/';
    self.pklPath = pklPath
    # three spectra namehead
    keHead = 'KE'
    divHead = 'DIV'
    zetaHead = 'ZETA'
    zetaGHead = 'ZETA_G'
    zetaAGHead = 'dZETA'
    lapPHead = 'Pk5'
    # make a list
    self.fieldHeadList = [keHead, divHead, zetaHead, zetaGHead, zetaAGHead, lapPHead]
    # three locations strings
    self.locStrList = ['surf','sekm', '-200m']

  def unzip_pickle(self, indLoc, indField, clmin, clmax):
    appendixStr = 'days.pkl'''
    specPKL = [['']*len(self.fieldHeadList)]*len(self.locStrList)
    fieldnameArrL = '_'.join([self.fieldHeadList[indField], \
                  self.locStrList[indLoc],self.casename,'{:d}'.format(int(self.tDuration)),appendixStr])
    specPKL = load_pickle(self.pklPath + fieldnameArrL,'rb')
    k1d = 0.5*(specPKL['binedge_kappa'][0:-1] + specPKL['binedge_kappa'][1:])*1000
    omg = specPKL['omg_eff']
    specMat = specPKL['KE_kappaT_eff']
    # return k1d,omg,specMat
    K,W=np.meshgrid(k1d,omg)
    plt.figure()
    plt.pcolormesh(k1d,omg,tch.log(tch.tensor(W)*tch.tensor(K)*specMat)/tch.log(tch.tensor(10.0)),cmap = 'terrain')
    plt.colorbar()
    plt.title(self.fieldHeadList[indField]+', '+self.locStrList[indLoc])
    plt.clim(clmin,clmax)
    plt.yscale('log'),plt.ylim([0.06,4])
    plt.xscale('log'),plt.xlim([0.004,k1d.max()/np.sqrt(2.0)])
    plt.suptitle(self.casename)

#%%
pathSpec = home_dir+'/postproc/img';

clmin = 4
clmax = 10

for indLoc in range(3):
  for indField in range(5):
    # del spectra_case
    spectra_case = spectraPKL('ex16_ucon_9', 36, pathSpec) # spectra_tp16uvar9.loadPKL();
    spectra_case.unzip_pickle(indLoc, indField, clmin, clmax)    
    
    del spectra_case
    spectra_case = spectraPKL('ex16_uvar_9', 36, pathSpec) # spectra_tp16uvar9.loadPKL()
    spectra_case.unzip_pickle(indLoc, indField, clmin, clmax)
    
    spectra_case = spectraPKL('tp16_uvar_9', 37, pathSpec)# spectra_tp16uvar9.loadPKL();
    spectra_case.unzip_pickle(indLoc, indField, clmin, clmax)
    
    del spectra_case
    spectra_case = spectraPKL('tp16_ucon_9', 36, pathSpec)# spectra_tp16uvar9.loadPKL();
    spectra_case.unzip_pickle(indLoc, indField, clmin, clmax)
