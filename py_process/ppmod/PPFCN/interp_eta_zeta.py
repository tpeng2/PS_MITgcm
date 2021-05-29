#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:33:41 2021

Interpolate 2D field between shifted coordinates

@author: tpeng2
"""

import numpy as np
from scipy import interpolate

def interp_eta_zeta(zeta,dx_move,dy_move,Lx,Ly,max_iters,min_res):
    Ny,Nx=zeta.shape
    dx=Lx/Nx
    dy=Ly/Ny
    x_orig=np.linspace(0,Lx-dx,num=Nx)
    y_orig=np.linspace(0,Ly-dy,num=Ny)
    # set interpolation framework based on input, quintic
    f_orig=interpolate.interp2d(x_orig,y_orig,zeta,kind='quintic')
    # now shifiting coordinates
    x_shift=np.linspace(dx_move,(Lx-dx)+dx_move,num=Nx)
    y_shift=np.linspace(dy_move,(Ly-dy)+dy_move,num=Ny)
    # first attemp of interpolation
    zeta_shift=f_orig(x_shift,y_shift)
    # set interpolation framework based on shifted coordinate
    f_shift=interpolate.interp2d(x_shift,y_shift,zeta_shift,kind='quintic')
    # shift back for checking the residual
    zeta_at_orig=f_shift(x_orig,y_orig)
    # Calculate difference
    dzeta_intp=zeta_at_orig-zeta
    # set residual
    residual = np.zeros(max_iters)
    residual[0] = np.max(dzeta_intp)
    print('starting residual = '+str(residual[0]))
    
    # Iterate to reduce residual (interpolate residual and then feedback)
    for i in range(max_iters):
        if (i>=1) and (residual[i-1] < min_res):
            print('Residual limit is satisfied. Last step: '+ str(residual[i-1])+'; Limit: '\
                  +str(min_res))
            return zeta_shift,residual;
        else:
            tmp_rand_df=np.random.randint(2, size=1)
            kind_df='linear';kind_shift='cubic';
            if np.mod(i,3)==0:
                if tmp_rand_df==1:
                    kind_df='linear'; 
                else:
                    kind_df='cubic';
            
            # define interpolation for difference
            df=interpolate.interp2d(x_orig,y_orig,dzeta_intp,kind=kind_df)
            # interpolate difference to the shifted coordinate
            dzeta_shift=df(x_shift,y_shift)
            # Set regression rate with a noise
            if i>1:
                if residual[i]>residual[i-1]:
                    damp=-0.5* residual[i]/residual[0]
                else:
                    damp=0.1* (residual[i]/residual[0])
            elif i==0:
                damp=0.05
            zeta_shift=zeta_shift-dzeta_shift*(1+np.random.randn(1)*damp)
            # redefine the shift-back interpolation
            f_shift=interpolate.interp2d(x_shift,y_shift,zeta_shift,kind=kind_shift)
            # shift back for checkin the residual
            zeta_at_orig=f_shift(x,y)
            dzeta_intp=zeta_at_orig-zeta
            residual[i]=np.max(np.abs(dzeta_intp))
            print('Step = '+str(i)+';','Res = '+str(residual[i]))
    return zeta_shift,residual,x_shift,y_shift