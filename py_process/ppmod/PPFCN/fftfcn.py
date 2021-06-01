#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 03:07:57 2021
PPFCNS.fftfcn
@author: tpeng
"""

import numpy as np
import torch as tch
import PPFCN.ppfcns as ppf


# class proc_2D_field:
#     """ To handle 2D fields, using FFT and other tools """
#     # Initialize datatype
#     def __init__(self,fpath,Nx,Ny):
#         """constructor"""
#         self.fpath=fpath
#         self.fhead=fhead
#         self.nx=Nx
#         self.ny=Ny
#         print('proc_2D_field: Information of path and filename is loaded.')
        
#     def read_2D_fields(self,):
    

""" Gaussian Filter (from internet) """ 
import math
import numbers
import torch
from torch import nn
from torch.nn import functional as F

class GaussianSmoothing(nn.Module):
    """
    Apply gaussian smoothing on a
    1d, 2d or 3d tensor. Filtering is performed seperately for each channel
    in the input using a depthwise convolution.
    Arguments:
        channels (int, sequence): Number of channels of the input tensors. Output will
            have this number of channels as well.
        kernel_size (int, sequence): Size of the gaussian kernel.
        sigma (float, sequence): Standard deviation of the gaussian kernel.
        dim (int, optional): The number of dimensions of the data.
            Default value is 2 (spatial).
    """
    def __init__(self, channels, kernel_size, sigma, dim=2):
        super(GaussianSmoothing, self).__init__()
        if isinstance(kernel_size, numbers.Number):
            kernel_size = [kernel_size] * dim
        if isinstance(sigma, numbers.Number):
            sigma = [sigma] * dim

        # The gaussian kernel is the product of the
        # gaussian function of each dimension.
        kernel = 1
        meshgrids = torch.meshgrid(
            [
                torch.arange(size, dtype=torch.float32)
                for size in kernel_size
            ]
        )
        for size, std, mgrid in zip(kernel_size, sigma, meshgrids):
            mean = (size - 1) / 2
            kernel *= 1 / (std * math.sqrt(2 * math.pi)) * \
                      torch.exp(-((mgrid - mean) / (2 * std)) ** 2)

        # Make sure sum of values in gaussian kernel equals 1.
        kernel = kernel / torch.sum(kernel)

        # Reshape to depthwise convolutional weight
        kernel = kernel.view(1, 1, *kernel.size())
        kernel = kernel.repeat(channels, *[1] * (kernel.dim() - 1))

        self.register_buffer('weight', kernel)
        self.groups = channels

        if dim == 1:
            self.conv = F.conv1d
        elif dim == 2:
            self.conv = F.conv2d
        elif dim == 3:
            self.conv = F.conv3d
        else:
            raise RuntimeError(
                'Only 1, 2 and 3 dimensions are supported. Received {}.'.format(dim)
            )

    def forward(self, input):
        """
        Apply gaussian filter to input.
        Arguments:
            input (torch.Tensor): Input to apply gaussian filter on.
        Returns:
            filtered (torch.Tensor): Filtered output.
        """
        return self.conv(input, weight=self.weight, groups=self.groups)




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

#%% Frequency filter
def fourier_freq_filter(A,Fs,f_cut_left,f_cut_right):
    Nfile=len(A)
    #set index
    if np.mod(Nfile//2,2)==0:
        nfreq=tch.linspace(0,Nfile//2-1,Nfile//2).type(tch.long)
    else:
        nfreq=tch.linspace(0,Nfile//2,Nfile//2+1).type(tch.long)
    Afilter=ppf.gen_filter_1d(A)
    Afilter=tch.tensor(Afilter[2])
    Afft=tch.fft.rfft(A*Afilter)[nfreq]
    spd=n_to_freq(gen_n_vec(Nfile),Fs)[nfreq]
    # remove unwanted signals
    nfreq_cut_left=tch.where(spd<f_cut_left)
    nfreq_cut_right=tch.where(spd>=f_cut_right)
    Afft[nfreq_cut_left]=tch.complex(tch.tensor(0.),tch.tensor(0.))
    Afft[nfreq_cut_right]=tch.complex(tch.tensor(0.),tch.tensor(0.))
    A_filtered=tch.fft.irfft(Afft,axis=0)
    return A_filtered.detach().clone()