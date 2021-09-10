#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 19:01:55 2020

@author: pmazumdar
"""

from astropy.io import fits
import numpy as np


#hdu_T = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_2.fits')[0]
hdu_T = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/fits_maps/12CO_fits/moment-maps/mom-0/G348.25_12CO_mom0.fits')[0]

data_T = hdu_T.data

#data_T[np.isnan(data_T)] = 0

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
#&%
#&% Removing the extra dimension
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&


def remove_hd_dim(hd, dim=3):

    """
    Keep information from datacube header
    only related to dimension to keep
    
    hd: astropy.header instance
    
    dim: int
    dimensions to keep, default=3
    
   """
    if dim == 2:
        nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS2'],hd['NAXIS1']])).header    
        dimlist = ['1','2']
    if dim == 3:
        nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS3'],hd['NAXIS2'],hd['NAXIS1']])).header
        dimlist = ['1','2','3']        
        
    for i in dimlist:
        for t in ['CRVAL','CRPIX','CDELT','CTYPE','CROTA','CUNIT']:
            if hd.get(t+i) != None:
                nhd[t+i] = hd[t+i]

    for t in ['BUNIT','BMAJ','BMIN','BPA','RESTFRQ']:
        if hd.get(t) != None:
            nhd[t] = hd[t]

    return nhd

data = data_T.squeeze()  # remove 4th dimension
hd = remove_hd_dim(hdu_T.header, dim=2) # change the header
hd['CTYPE1'] = 'GLON-GLS'
hd['CTYPE2'] = 'GLAT-GLS'
hdu_T.data = data
hdu_T.header = hd
hdu_T.writeto('/home/pmazumdar/Documents/LASMA/Reduction/fits_maps/12CO_fits/moment-maps/mom-0/G348.25_12CO_mom0.fits', overwrite=1)
