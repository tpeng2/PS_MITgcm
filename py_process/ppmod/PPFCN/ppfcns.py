#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 23:24:35 2021

Post Process functions and classes

@author: tpeng
"""





import numpy as np
import torch as tch
import glob
import matplotlib.pyplot as plt
import MITgcmutils as utils
import math
#%% Model parameters
def gen_f_coriolis(f_0,beta,x_vec,y_vec):
    f=f_0+beta*y_vec
    f2D=np.tile(f,(len(x_vec),1)).transpose()
    return f,f2D



# read vertical grids

def load_vertical_grids(casepath):
    if casepath[-1]!='/':
        casepath+='/'
    RC=utils.rdmds(casepath+'RC')
    RF=utils.rdmds(casepath+'RF')
    DRC=utils.rdmds(casepath+'DRC')
    DRF=utils.rdmds(casepath+'DRF')
    return RC,RF,DRC,DRF
def load_horizontal_grids(casepath):
    if casepath[-1]!='/':
        casepath+='/'
    XC=utils.rdmds(casepath+'XC') 
    YC=utils.rdmds(casepath+'YC') 
    XG=utils.rdmds(casepath+'XG') 
    YG=utils.rdmds(casepath+'YG')
    return XC,YC,XG,YG        
#%% File process
# Search file name
def search_mdsdata(fpath,fhead):
    fname=sorted(glob.glob(fpath+fhead+'*data',recursive=True));
    return fname

def search_binary(fpath,fhead):
    fname=sorted(glob.glob(fpath+fhead,recursive=True));
    return fname

# Trim File name given 

def trim_fname(filename,suffix,start_id,end_id):
    suffix_len=len(suffix)
    file_head=filename[0][:-(10+suffix_len)]
    nfiles=len(filename)
    file_index=np.zeros(nfiles)
    for i in np.arange(nfiles):
        if suffix_len == 0:
            file_index[i]=filename[i][-(10+suffix_len):]
        else:
            file_index[i]=filename[i][-15:-5];
    ind_remain=np.intersect1d(np.where(file_index<=end_id),np.where(file_index>=start_id))
    return filename[slice(ind_remain[0],ind_remain[-1]+1,1)]

# get 1D FFT filter
def gen_filter_1d(A):
    dim_A=len(A.shape)
    if(dim_A>3):
        print('error in dimensions!')
    NA=A.shape
    flter_tmp=np.zeros(dim_A,dtype=np.object)
    for i in np.arange(dim_A):
        axis_array=np.arange(NA[i])
        flter_tmp[i]=0.5*(1-np.cos(2*np.pi/NA[i]*axis_array))
    return flter_tmp

# Read binary 2D slides
def read_2Dfield(fpath,fhead,fid,Ny,Nx):
    fid_str=format(fid,'010d')
    M_read=np.fromfile(fpath+'/'+fhead+'.'+fid_str,dtype=np.float32).reshape(Ny,Nx)
    return M_read

#% Consolidate 2D fields to a single 3D matirx
def consolidate_2Dfield(fpath_2D,fheader,start_ind,end_ind,Ny,Nx,Nfiles):
    print('start reading ' + fheader + ' fields')
    M_2D_fname = search_binary(fpath_2D,fheader+'*')
    M_2D_fname_trimmed = trim_fname(M_2D_fname,'',start_ind,end_ind)
    Nfiles=len(M_2D_fname_trimmed);
    file_ind=[]
    M_2Dt=np.zeros((Nfiles,Ny,Nx))
    for i in range(Nfiles):
        file_ind.append(M_2D_fname_trimmed[i][-10:])
        M_2Dt[i,:,:]=np.fromfile(fpath_2D+fheader+'.'+file_ind[i],dtype=np.float32).reshape(Ny,Nx)
    return M_2Dt


# plot 2D functions
    
def plot_2d_colorbar(imgfiled,xvec,yvec,xlabel_str,ylabel_str,title_str,cmap_str,clim_range,plt_asp=1,fname=None,ftype=None,fpath=None):
    from datetime import datetime
    fig = plt.figure();
    # ax = fig.add_subplot(111,aspect=len(yvec)/len(xvec));
    ax = fig.add_subplot(111,aspect=plt_asp);
    plm=ax.pcolormesh(xvec,yvec,imgfiled,cmap=plt.cm.get_cmap(cmap_str),shading='auto');
    fig.colorbar(plm);plm.set_clim(clim_range);
    ax.set_xlabel(xlabel_str,fontsize='small');
    ax.set_ylabel(ylabel_str,fontsize='small');
    ax.set_title(title_str,fontsize='small')
    # ax.axis('tight')
    if fname==None:
        return
    elif ftype==None and fpath == None:
        ftype='.png'
        fpath='./'
    fig.set_size_inches(8,6)
    fig.savefig(fpath+'/'+fname+ftype, dpi=300, bbox_inches='tight',pad_inches = 0)
    del fig


#%% Make U, V vorticity
# dim [Ny,Nx]
def get_vorticity(U,V,dx,dy):
    Nx=U.shape[2]; Ny=V.shape[1];
    dUsurf_dy=np.diff(U,n=1,axis=1)/dy # ==> (Ny-1,Nx) at zeta point
    dVsurf_dx=np.diff(V,n=1,axis=2)/dx # ==> (Ny, Nx-1) at zeta point
    # dUsurf_dx=np.diff(U,n=1,axis=0)/dx # ==> (Ny-1,Nx) at zeta point
    # dVsurf_dy=np.diff(V,n=1,axis=1)/dy # ==> (Ny, Nx-1) at zeta point
    # Add edge values:
    if len(U.shape) == 2:
        # 2D
        dU_dy_ff=np.zeros([Ny,Nx])
        dV_dx_ff=np.zeros([Ny,Nx])
        dU_dx_ff=np.zeros([Ny,Nx])
        dV_dy_ff=np.zeros([Ny,Nx])
    elif len(U.shape)==3:
        Nfile_UV=U.shape[0]
        # Add edge values:
        dU_dy_ff=np.zeros([Nfile_UV,Ny,Nx])
        dV_dx_ff=np.zeros([Nfile_UV,Ny,Nx])
        dU_dx_ff=np.zeros([Nfile_UV,Ny,Nx])
        dV_dy_ff=np.zeros([Nfile_UV,Ny,Nx])
    dU_dy_ff[:,1:,:]=dUsurf_dy
    dV_dx_ff[:,:,1:]=dVsurf_dx
    # dU_dx_ff[1:,:]=dUsurf_dx
    # dV_dy_ff[:,1:]=dVsurf_dy
    
    # edge value
    dU_dy_ff[:,0,:]=U[:,0,:]/(dy/2)
    dV_dx_ff[:,:,0]=(V[:,:,0]-V[:,:,-1])/(dx)
    # dV_dy_ff[0,:]=V[0,:]/(dy/2)
    # dU_dx_ff[:,0]=(U[:,0]-U[:,-1])/(dx)
    # vorticity
    zeta=dV_dx_ff-dU_dy_ff;
    # div=dV_dy_ff+dU_dx_ff
    return zeta#,div




#%% Laplacian of 2D matrix

def calc_Laplacian(M,dx,dy):
    dM=np.gradient(M);
    dM_dx=dM[1]/dx; ddM_dx2=np.gradient(dM_dx)[1]/dx;
    dM_dy=dM[0]/dy; ddM_dy2=np.gradient(dM_dy)[0]/dy;
    StrM=ddM_dx2-ddM_dy2;
    LapM=ddM_dx2+ddM_dy2;
    return LapM,StrM

#%% Interpolate to shifted coordinates
    
#                   >>> dV/dx first diff value
#     |       |      |       |  #%%
#    zeta --  v  -- zeta --  v  -- zeta --
#     |       |      |       |       |    
#     u      eta     u      eta      u    
#     |       |      |       |       |    
#    zeta --  v  -- zeta --  v  -- zeta --
#     |       |      |       |       |    
#     u      eta     u      eta      u    
#     |       |      |       |       |    
#    zeta --  v  -- zeta --  v  -- zeta -- <<< dU/dy first diff value
#     |       |      |       |       |    
#     u      eta     u      eta      u   
#     |       |      |       |       |    
#    zeta --  v  -- zeta --  v  -- zeta --

# Overlapped field between dU/dy and dV/dx (indices start with 0):
# dU/dy ==> [:,1:Nx-1]
# dV/dx ==> [1:Ny-1,:]

from scipy import interpolate

def interp_eta_zeta(zeta,dx_move,dy_move,Lx,Ly,max_iters,min_res_rto):
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
        if (i>=1) and (residual[0]/residual[i-1]<min_res_rto):
            print('Residual limit is satisfied. Last step: '+ str(residual[i-1])+'; Limit: '\
                  +str(min_res_rto))
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
            zeta_at_orig=f_shift(x_orig,y_orig)
            dzeta_intp=zeta_at_orig-zeta
            residual[i]=np.max(np.abs(dzeta_intp))
            print('Step = '+str(i)+';','Max. Res. = '+str(residual[i]))
    if dx_move>0:
        x_output=x_shift
    else:
        x_output=x_orig
    
    if dy_move>0:
        y_output=y_shift
    else:
        y_output=y_orig
            
    return zeta_shift,residual,x_output,y_output

#%%
T_ref_initial=[7.950849079316678,
   7.861182766074640,
   7.758501159111756,
   7.644014320341626,
   7.518949428341353,
   7.384534361378004,
   7.241982877054776,
   7.092481545507864,
   6.937178485393330,
   6.777173921599222,
   6.613512496273191,
   6.447177291937836,
   6.279085460547976,
   6.110085285318365,
   5.940954654392985,
   5.772400601665104,
   5.605059946543879,
   5.439500785227247,
   5.276224719089327,
   5.115669619011567,
   4.958212902765920,
   4.804175099171582,
   4.653823654602955,
   4.507376872268201,
   4.365007914840069,
   4.226848754558513,
   4.092994087248265,
   3.963505123811641,
   3.838413204206522,
   3.717723231270353,
   3.601416939778244,
   3.489455892813602,
   3.381784295548079,
   3.278331575049326,
   3.179014745505644,
   3.083740540667922,
   2.992407345986870,
   2.904906924762290,
   2.821125982040719,
   2.740947511325965,
   2.664252004407315,
   2.590918502271465,
   2.520825511848193,
   2.453851784398542,
   2.389876984812668,
   2.328782256966866,
   2.270450691464814,
   2.214767728127219,
   2.161621461281368,
   2.110902882070932,
   2.062506091514641,
   2.016328423160970,
   1.972270554467213,
   1.930236564949263,
   1.890133966363022,
   1.851873705679003,
   1.815370148096710,
   1.780541031737083,
   1.747307418892613,
   1.715593625274912,
   1.685327139696111,
   1.656438546213544,
   1.628861424414074,
   1.602532256212109,
   1.577390327087151,
   1.553377627905139,
   1.530438749594961,
   1.508520789119703,
   1.487573241811485,
   1.467547912093703,
   1.448398810469260,
   1.430082068765359,
   1.412555845156073,
   1.395780241814092,
   1.379717220714235,
   1.364330523524627,
   1.349585601652287,
   1.335449539276355,
   1.321890992789812,
   1.308880123744298,
   1.296388542007769,
   1.284389252821320,
   1.272856605327441,
   1.261766246736765,
   1.251095079535491,
   1.240821227616209,
   1.230923999452373,
   1.221383860380165,
   1.212182408931796,
   1.203302355299647,
   1.194727504350672,
   1.186442745994009,
   1.178434046958083,
   1.170688449023884,
   1.163194068401473,
   1.155940105994612,
   1.148916857212706,
   1.142115727922288,
   1.135529257917519,
   1.129151147416661,]