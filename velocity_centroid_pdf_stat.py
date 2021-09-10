#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 13:08:19 2020

@author: pmazumdar

Task: Calculate the velocity-centroid pdf from the moment-1 map

"""
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import hist
import numpy as np
import matplotlib
#from astropy.modeling import models, fitting
#from scipy.optimize import curve_fit


hdu_mom1 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_1.fits')[0]
snr = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_snr.fits')[0]

data = hdu_mom1.data
#data = data[snr.data>5]


flat_data = data.flatten()
flat_snr = snr.data.flatten()

#flat_data[flat_data < -60] = np.nan
flat_data[flat_data > -10] = np.nan

flat_snr = flat_snr[~np.isnan(flat_data)]
flat_data = flat_data[~np.isnan(flat_data)]

flat_data = flat_data[flat_snr>5]
flat_snr = flat_snr[flat_snr>5]

flat_snr[flat_snr>100] = 100

matplotlib.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1, 1, figsize=(6.5, 6))
fig.subplots_adjust(bottom=0.15)

my_hist = hist(flat_data, bins=100, weights=flat_snr, histtype='step',ax=ax,density=True,color='k')

bin_centers = np.array([])
my_data = np.array(my_hist[0])

for i in np.arange(len(my_hist[1])-1):
    bin_centers = np.append(bin_centers,(my_hist[1][i+1]+my_hist[1][i])/2.0)

bin_width = bin_centers[1]-bin_centers[0]
nevents = len(flat_data)

#center_error = np.sqrt(my_data.max())
#rel_error = (np.abs(my_bin-bin_center)/10)**2
#my_error = center_error*my_data*rel_error
my_data = my_data*nevents*bin_width
my_error = np.sqrt(my_data)
#my_error[my_error<2] = 0.5
my_error = my_error/nevents/bin_width
#my_error = 2*my_error
my_data = my_data/nevents/bin_width

ax.errorbar(bin_centers,my_data,yerr=my_error,ecolor='k',fmt='k',linestyle = 'None',linewidth=1.2)
ax.set_yscale('log')
ax.set_xlabel('centroid velocity [km.s$^{-1}$]')
ax.set_ylabel(r'P(v) [(km.s$^{-1}$)$^{-1}$]')

plt.tight_layout()
plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/v_centroid_pdf.eps',dvi=200)

v = bin_centers
P = my_data

dv = v[2]-v[1]
v_bar = np.sum(dv*P*v)
var_v = np.sum(dv*P*((v-v_bar)**2))
sig_v = np.sqrt(var_v)
kurt_v = np.sum(dv*P*((v-v_bar)**4))/(var_v**2)

errors = (dv/v)+(my_error/P)
errors[np.isnan(errors)]=0.0
error_v_bar = np.sqrt(np.sum(errors**2))

errors = 2*((np.sqrt(((error_v_bar/v_bar)**2+(dv/v)**2)))/(v-v_bar))+(my_error/P)
errors[np.isnan(errors)]=0.0
error_var_v = np.sqrt(np.sum(errors**2))

error_sig_v = 0.5*sig_v*(error_var_v/var_v)

errors = 4*((np.sqrt(((error_v_bar/v_bar)**2+(dv/v)**2)))/(v-v_bar))+(my_error/P)
errors[np.isnan(errors)]=0.0
error_k = np.sqrt(np.sum(errors**2))
error_kurt_v = (2*error_var_v/var_v)+(error_k/np.sum(dv*P*((v-v_bar)**4)))

print(r'Mean = {:.1f}+/-{:.1f} , std = {:.1f}+/-{:.1f}, K = {:.1f}+/-{:.1f}'.\
      format(v_bar,error_v_bar,sig_v,error_sig_v,kurt_v,error_kurt_v))