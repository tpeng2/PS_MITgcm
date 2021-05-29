# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 15:27:56 2021

Generate New initial case files

@author: tpeng
"""

import numpy as np
import torch as tch
import matplotlib.pyplot as plt
import numpy.matlib as matlib 
def gen_steady_windstress(dx,dy,Lx,Ly,tau_0):
    Nx=int(Lx/dx)
    Ny=int(Ly/dy)
    xvec=np.linspace(dx/2,Lx-dx/2,num=int(Nx))
    yvec=np.linspace(dy/2,Ly-dx/2,num=int(Ny))
    tau=tau_0*(np.sin(np.pi*np.arange(len(yvec))/Ny))**2
    # make it 2D, homogeneous in x-direction
    tau_2Dxy=matlib.repmat(tau,Nx,1).squeeze()
    plt.pcolormesh(yvec/1e3,xvec/1e3,tau_2Dxy,shading='auto'),plt.colorbar(),plt.xlabel('y[km]');plt.ylabel('x[km]')
    plt.title('wind stress ([Nx,Ny])')
    return tau_2Dxy

#%% Test bed
