#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:37:47 2021

@author: tpeng
"""

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

args = parser.parse_args()

path_tsr=args.path_tsr
fname_vort_tsr=args.fname_vort_tsr
fname_vort_eta_tsr=args.fname_vort_eta_tsr
Fs=int(args.Fs)
dt_sec=int(args.dt_sec)

#%%

path_tsr=home_dir+'/postproc/tsr/'
fname_vort_tsr='fltd_high_Vort_surf_rel_uvar_9_0000480118_0000549958.pt'
fname_vort_eta_tsr='fltd_high_Vort_eta_surf_rel_uvar_9_0000480118_0000549958.pt'

vort=tch.load(path_tsr+fname_vort_tsr)
vort_eta=tch.load(path_tsr+fname_vort_eta_tsr)
#%%

dx=500
dy=500
Lx=480e3
Ly=960e3
xc=np.linspace(dx/2,Lx-dx/2,int(Lx//dx))
yc=np.linspace(dy/2,Ly-dy/2,int(Ly//dy))

#%%
