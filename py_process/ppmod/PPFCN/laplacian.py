#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:46:10 2021

Calculate Laplacian
https://discuss.pytorch.org/t/how-to-calculate-laplacian-sum-of-2nd-derivatives-in-one-step/41667/3

@author: tpeng2
"""

import torch
import numpy as np

def calc_Laplacian(M,dx,dy):
    dM=np.gradient(M);
    dM_dx=dM[1]/dx; ddM_dx2=np.gradient(dM_dx)[1]/dx;
    dM_dy=dM[0]/dy; ddM_dy2=np.gradient(dM_dy)[0]/dy;
    StrM=ddM_dx2-ddM_dy2;
    LapM=ddM_dx2+ddM_dy2;
    return LapM,StrM