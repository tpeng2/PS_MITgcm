#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 03:11:24 2021
__init__.py
@author: tpeng
"""


from scipy import signal
import numpy as np
import torch as tch
import torch.fft as tchfft
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import MITgcmutils as utils
#import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import scipy.io as scipyio
import glob
import PPFCN.ppfcns as ppf


#%% Class
class cls_Ext_xy_2Dmap:
    # Initialization
    def __init__(self,dx,dy,Lx,Ly,fhead,ext_savepath,groupname,casename,start_ind,end_ind,ftype=''):
        self.dx=dx
        self.dy=dy
        self.Lx=Lx
        self.Ly=Ly
        self.ny=int(Ly/dy)
        self.nx=int(Lx/dx)  
        self.xc=np.linspace(dx/2,Lx-dx/2,self.nx)
        self.yc=np.linspace(dy/2,Ly-dy/2,self.ny)
        self.fhead=fhead
        self.__ext_savepath=ext_savepath
        self.__groupname=groupname
        self.__casename=casename
        self.__start_ind=start_ind
        self.__end_ind=end_ind
        self.__ftype=ftype
    # Assign path and file index range
    def assign_path(self):
        if self.__ext_savepath[-1]!='/':
            self.__ext_savepath=self.__ext_savepath+'/'
        self.path_ext=self.__ext_savepath+self.__groupname+'/'+self.__casename+'/'
        rawfname=ppf.search_binary(self.path_ext,self.fhead+'*')
        # trim file into desired range
        self.fname=ppf.trim_fname(rawfname,self.__ftype,self.__start_ind,self.__end_ind)
        print('file name are trimmed now')
        self.Nfile=len(self.fname)
        # Initialize 2D
        self.mat=np.zeros((self.Nfile,self.ny,self.nx))
        del rawfname
    # Load files
    def load_ext_files(self):
        file_ind=[]
        for i in range(self.Nfile):
            if len(self.__ftype)==0:
                file_ind.append(self.fname[i][-(10):])
            else:
                file_ind.append(self.fname[i][-(10+len(self.__ftype)):-len(self.__ftype)])
            self.mat[i,:,:]=np.fromfile(self.fname[i],dtype=np.float32).reshape(self.ny,self.nx) 
        self.file_ind=file_ind
        del file_ind

#%% Secondary class

class Exf_UVP_2D:
    def __init__(self,dx,dy,Lx,Ly,U_fhead,V_fhead,P_fhead,ext_savepath,groupname,casename,start_ind,end_ind):
        self.__dx=dx
        self.__dy=dy
        self.__Lx=Lx
        self.__Ly=Ly
        self.__ext_savepath=ext_savepath
        self.__Utmp=cls_Ext_xy_2Dmap(self.__dx,self.__dy,self.__Lx,self.__Ly,U_fhead,self.__ext_savepath,groupname,casename,start_ind,end_ind)
        self.__Vtmp=cls_Ext_xy_2Dmap(self.__dx,self.__dy,self.__Lx,self.__Ly,V_fhead,self.__ext_savepath,groupname,casename,start_ind,end_ind)
        self.__Ptmp=cls_Ext_xy_2Dmap(self.__dx,self.__dy,self.__Lx,self.__Ly,P_fhead,self.__ext_savepath,groupname,casename,start_ind,end_ind)
        self.__xc=self.__Utmp.xc
        self.__yc=self.__Utmp.yc
        self.f,self.f_2D=ppf.gen_f_coriolis(-1e-4, 1.5e-11, self.__xc, self.__yc)
        self.__nx=self.__Lx/self.__dx
        self.__ny=self.__Ly/self.__dy
        #

    def load_P(self):
        self.__Ptmp.assign_path()
        self.__Ptmp.load_ext_files()
        self.P=self.__Ptmp.mat
        self.P_filters=ppf.gen_filter_1d(self.P)
        self.P_fileind=self.__Ptmp.file_ind 
        
    def calc_vort_eta(self,niter_max):
        self.Vort_eta=tch.zeros(self.P.shape)
        # self.Str_P_at_zeta=np.zeros(self.P.shape)
        for i in range(len(self.P)):
            print('Interpolate file #'+str(i)+' of '+str(len(self.P)),end='')
            eta_at_zeta,res_eta_max,x_zeta,y_zeta=ppf.interp_eta_zeta(self.P[i]/9.81,-self.__dx//2,-self.__dy//2,self.__Lx,self.__Ly,niter_max,2e-6)
            f,f_2D=ppf.gen_f_coriolis(-1e-4, 1.5e-11, x_zeta, y_zeta)
            Lap_P_at_zeta,Str_P_at_zeta=ppf.calc_Laplacian(eta_at_zeta, self.__dx, self.__dy)
            self.Vort_eta[i,:,:]=tch.tensor(Lap_P_at_zeta/(f_2D*f_2D)*9.81)
            
    def calc_geostrophic_velocity(self):
            self.ugeo=-np.gradient(self.P,axis=1)/self.__dy/(self.f_2D) # P is P/rho
            self.vgeo=np.gradient(self.P,axis=2)/self.__dx/(self.f_2D) # P is P/rho
            
    def load_U(self):
        self.__Utmp.assign_path()
        self.__Utmp.load_ext_files();
        self.U=tch.tensor(self.__Utmp.mat)
        self.U_filters=ppf.gen_filter_1d(self.U)
        self.U_fileind=self.__Utmp.file_ind 
    def load_V(self):
        self.__Vtmp.assign_path()
        self.__Vtmp.load_ext_files()
        self.V=tch.tensor(self.__Vtmp.mat)
        self.V_filters=ppf.gen_filter_1d(self.V)
        self.V_fileind=self.__Vtmp.file_ind 

    def calc_rel_vorticity(self):
        self.vort=tch.tensor(ppf.get_vorticity(self.U,self.V,self.__dx,self.__dy)).type(tch.float)
    
    def calc_rotary_spectrum(self):
        self.U_filters=ppf.gen_filter_1d(self.U)
        self.V_filters=ppf.gen_filter_1d(self.V)
        U_tensor=tch.tensor(self.U_filters[0][:,None,None])*self.U
        V_tensor=tch.tensor(self.V_filters[0][:,None,None])*self.V
        R=tch.complex(U_tensor.type(tch.float),V_tensor.type(tch.float))
        del U_tensor; del V_tensor;
        fft_R=tch.fft.fft(R,dim=0)
        del R
        self.PRR=(fft_R)*tch.conj(fft_R)
        del fft_R
        
    def calc_rotary_spectrum_geo_residual(self):
        U_tensor=tch.tensor(self.U_filters[0][:,None,None])*(self.U-self.ugeo)
        V_tensor=tch.tensor(self.V_filters[0][:,None,None])*(self.V-self.vgeo)
        R=tch.complex(U_tensor.type(tch.float),V_tensor.type(tch.float))
        del U_tensor; del V_tensor;
        fft_R=tch.fft.fft(R,dim=0)
        del R
        self.PRR_dgeo=(fft_R)*tch.conj(fft_R)
        del fft_R
        
    def calc_rotary_spectrum_geo(self):
        self.calc_geostrophic_velocity()
        self.U_filters=ppf.gen_filter_1d(self.ugeo)
        self.V_filters=ppf.gen_filter_1d(self.vgeo)
        U_tensor=tch.tensor(self.U_filters[0][:,None,None])*(self.ugeo)
        V_tensor=tch.tensor(self.V_filters[0][:,None,None])*(self.vgeo)
        R=tch.complex(U_tensor.type(tch.float),V_tensor.type(tch.float))
        del U_tensor; del V_tensor;
        fft_R=tch.fft.fft(R,dim=0)
        del R
        self.PRR_geo=(fft_R)*tch.conj(fft_R)
        del fft_R
    

    def gen_ft_array(self,res_smp,len_smp):
        if np.mod(len_smp,2)==0:
            ftarr=np.linspace(0,(len_smp)/2/(res_smp*len_smp),int(len_smp)//2+1) #f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
        else:
            ftarr=np.linspace(0,(len_smp-1)/2/(res_smp*len_smp),int(len_smp)//2+1) #f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
        return ftarr
    
    
    def freq_filtering(self,dt_smp,freq_low,freq_high):
        self.omg=self.gen_ft_array(dt_smp,self.__Utmp.Nfile)
        omg_tensor=tch.tensor(self.omg)
        omg_drop_low=tch.where(omg_tensor<freq_low)
        omg_drop_high=tch.where(omg_tensor>=freq_high)
        if hasattr(self,'vort'):
            print('Filtering relative vorticity from '+str(freq_low)+' to '+str(freq_high))
            vort_fft=tch.fft.rfft(tch.tensor(self.U_filters[0][:,None,None])*self.vort,dim=0)
            vort_fft[omg_drop_low]=0
            vort_fft[omg_drop_high]=0
            self.vort_filtered=tch.fft.irfft(vort_fft,dim=0)
            if hasattr(self,'Vort_eta'):
                print('Filtering Laplacian SSH from '+str(freq_low)+' to '+str(freq_high))
                vort_fft=tch.fft.rfft(tch.tensor(self.U_filters[0][:,None,None])*self.Vort_eta,dim=0)
                vort_fft[omg_drop_low]=0
                vort_fft[omg_drop_high]=0
                self.vort_eta_filtered=tch.fft.irfft(vort_fft,dim=0)
        else:
            print('Relative vorticity has not been calculated yet.')
        
        
        
    def del_U(self):
        del self.U
    def del_V(self):
        del self.V
    def del_P(self):
        del self.P        
