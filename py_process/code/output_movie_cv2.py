#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 21:16:39 2021

@author: tpeng
"""

import cv2
import numpy as np
import glob
import sys
import os

home_dir=os.path.expanduser("~")
img_array = []
for filename in glob.glob(home_dir+'/postproc/img/moviefiles/freq_sekm_filter_vort_real_uvar_9_60_day_fr*'):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)
 
 
out = cv2.VideoWriter('freq_filter_vort_real_uvar_9_60_day.avi',cv2.VideoWriter_fourcc(*'DIVX'), 12, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()