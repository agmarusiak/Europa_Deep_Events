#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:49:10 2021

@author: marusiak
"""



from obspy import read, read_inventory
from obspy.io.xseed import Parser
from obspy.signal import PPSD
import numpy
import scipy
import matplotlib.pyplot as plt
import os
import glob

from obspy.imaging.cm import pqlx


fileZ="/Users/marusiak/Documents/GitHub/EuropaNoise/noise_records/ice20.pref_cat1.MXZ"
st=read(fileZ)
tr=st.select(channel='MXZ')[0]
print(tr.stats)

paz = {'gain': 1.0,
       'poles': [0j],
       'zeros': [],
       'sensitivity': 1.0}



ppsd=PPSD(tr.stats,paz,db_bins=(-325, -100, 1.0),)
ppsd.add(st)
median=ppsd.get_percentile(percentile=50)
print(median)
fname='/Users/marusiak/Documents/GitHub/EuropaNoise/noise_records/ice5_pref_cat20_MXZ.mat'
import scipy.io as sio
#fname4="PDFS/"+station+"data.mat"
sio.savemat(fname,{'median':median})
ppsd.plot(cmap=pqlx,period_lim=(0.3,600))