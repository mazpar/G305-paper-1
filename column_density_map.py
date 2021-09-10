#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:40:24 2020

@author: pmazumdar

Task: create a column density cube given optical depths and excitation temperature
"""
from astropy.io import fits
from astropy import wcs
import numpy as np
import aplpy

#T_map = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/temp_12.fits')[0]
T_map = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/temp.fits')[0]
tau_map = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/tau.fits')[0]

T_data = T_map.data
T_data[np.isnan(T_data)] = 0.0
T_data[T_data<0] = 0.0

hd = T_map.header
w = wcs.WCS(T_map.header)

tau_data = tau_map.data
tau_data[np.isnan(tau_data)] = 0.0
tau_data[tau_data<0] = 0.0


h = 6.62607004e-34
nu = 330.5880e9
k_B = 1.38064852e-23
factor = 1.1037108322258952e+25 #(8*pi/c**3)*(g2/g3)*(nu**3/A_32); velocity taken in km/s
k_hb = 0.3550599734823035 # k_B/h*B value in K^-1
vres = 0.5


n13_data = 1e-10*((factor*(1/(1-np.exp((-h*nu)/(k_B*T_data))))*(tau_data*abs(vres))))

z = k_hb*(T_data+(1/(3*k_hb)))  # calculate the partition function for each voxel

ntot_data = n13_data*(z/5)*np.exp(6/(T_data*k_hb))

ntot_data[np.isnan(ntot_data)] = 0.0
ntot_data[ntot_data<0.0] = 0.0

nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS3'],hd['NAXIS2'],hd['NAXIS1']])).header
dimlist = ['1','2','3']

for i in dimlist:
    for t in ['CRVAL','CRPIX','CDELT','CTYPE','CROTA','CUNIT']:
        if hd.get(t+i) != None:
            nhd[t+i] = hd[t+i]

for t in ['BMAJ','BMIN','BPA','RESTFRQ']:
    if hd.get(t) != None:
        nhd[t] = hd[t]

nhd['BUNIT'] = 'Jy/beam'


fits.writeto('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal.fits',ntot_data,header=nhd,overwrite=True)
