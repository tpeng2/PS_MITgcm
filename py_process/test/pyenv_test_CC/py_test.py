# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:15:38 2021
Make wave-speactra for 
@author: tpeng2
"""


import numpy as np
import torch as tch
import torch.fft as tchfft
import matplotlib
import matplotlib.pyplot as plt
import MITgcmutils as utils
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import scipy.io as scipyio


#%% define functions
output_path='./'
tch.pi=tch.tensor(np.pi)
pi=tch.pi
print(pi)
np.savetxt(output_path+'test_results.csv',(pi,pi,pi),delimiter=',')
