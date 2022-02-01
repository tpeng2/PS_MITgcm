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
    zetaGHead = 'ZETA_geo'
    zetaAGHead = 'ZETA_ageo'
    lapPHead = 'Lap_P'
    # make a list
    self.fieldHeadList = [keHead, divHead, zetaHead, zetaGHead, zetaAGHead, lapPHead]
    # three locations strings
    self.locStrList = ['surf','sekm', '-200m']

  def loadPKL(self):
    # construct file name
    appendixStr = 'days.pkl'
    self.fieldnameArrL = [['']*len(self.fieldHeadList)]*len(self.locStrList)
    self.spectra = [['']*len(self.fieldHeadList)]*len(self.locStrList)
    # format: field_caseame_appendix
    ifield = 0
    for field in self.fieldHeadList:
      iloc = 0
      for loc in self.locStrList:
        self.fieldnameArrL[iloc][ifield] = '_'.join([self.fieldHeadList[ifield], \
                                                     self.locStrList[iloc],self.casename,'{:d}'.format(int(self.tDuration)),appendixStr])
        self.specPKL[iloc][ifield] = load_pickle(self.pklPath + self.fieldnameArrL[iloc][ifield],'rb')
        iloc+=1;
      ifield+=1;

#%%
pathSpec = home_dir+'/postproc/img';

spectra_tp16uvar9 = spectraPKL('tp16_uvar_9',37,pathSpec)
spectra_tp16uvar9.loadPKL();