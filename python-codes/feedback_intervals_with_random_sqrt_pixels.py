#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:17:58 2020

@author: pmazumdar

Task: read Paper-1 section 7.3
"""

################# Plotting Average Spectra over a masked cube #################
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
import matplotlib
import aplpy
import scipy.stats as ss

################################
## Want to save output plots? ##
################################
save = False




############################
## Read the spectral cube ##
############################

cube = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_resample.fits')
#cube = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO_resample.fits')

nii_cube = cube.with_spectral_unit(u.km/u.s)
nii_subcube = nii_cube.spectral_slab(-60*u.km/u.s,0*u.km/u.s)




############################
##    Read the 8um Map    ##
############################

hdu_glm8 = fits.open('/home/pmazumdar/Documents/LASMA/Ancillary_Data/GLIMPSE/0.6_mosaics_v2.0/'+\
                     'GLM_30550+0000_mosaic_I4.fits')[0]
data_glm8 = hdu_glm8.data
w_glm8 = wcs.WCS(hdu_glm8.header)




############################
##  Read the moment Map   ##
############################
hdu_12CO = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_0.fits')[0]
hdu_mom1 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_1.fits')[0]




#########################################
##          Reprojections              ##
#########################################

from reproject import reproject_interp

#hdu_T = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal_sum_squeeze.fits')[0]
hdu_T = nii_subcube[0,:,:].hdu


array, footprint = reproject_interp(hdu_glm8, hdu_T.header)
array0, footprint0 = reproject_interp(hdu_12CO, hdu_T.header)  # reprojects moment maps on data cube for same dimensions 
array1,footprint1 = reproject_interp(hdu_mom1, hdu_T.header)


fits.writeto('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits', array,\
             hdu_T.header, clobber=True)
fits.writeto('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO-moment0.fits', array0,\
             hdu_T.header, clobber=True)
fits.writeto('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO-mom1.fits', array1,\
             hdu_T.header, clobber=True)




#########################################
##         Load Reprojections          ##
#########################################

hdu_dust = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
w_dust = wcs.WCS(hdu_dust.header)

hdu_12CO = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO-moment0.fits')[0]
w_12CO = wcs.WCS(hdu_12CO.header)

hdu_mom1 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO-mom1.fits')[0]




# ------> 12CO intensity based mask to avoid UC H II regions <----------

mask_nan = ~np.isnan(hdu_12CO.data) # only include non nan pixels
#mask_uchii = hdu_12CO.data<100    # upper limit on intensity to remove uc hii regions
mask_12CO = mask_nan


nii_subcube = nii_subcube.with_mask(mask_12CO)   # apply mask to spectral cube




#################################################################
# Aligning spectra based on moment 1 map / from velocity peaks ##
#################################################################

#----> mom-1 based shifts <----
mom1_data = hdu_mom1.data
# mom1_data[np.isnan(mom1_data)] = -40
mom1_data[mom1_data<-50] = -40 # values <-60 are pulled back to -45 since they are most likely from noise
mom1_data[mom1_data>-15] = -40 # values >-10 are pulled back to -45 since they are most likely from noise

mom1_shift = 2.0*np.ceil(-40.0-mom1_data) # Determining the channel number shift to use in numpy.roll
mom1_shift[np.isnan(mom1_shift)] = 0.0

peak_shift = mom1_shift.astype(int)
peak_shift[peak_shift < 0] += nii_subcube.shape[0]  # Use always a -ve shift, s.t. column_indices are valid.

i, j, k_indices = np.ogrid[:nii_subcube.shape[0], :nii_subcube.shape[1], :nii_subcube.shape[2]]  #re-indexing
i = i - peak_shift[np.newaxis,:,:]

aligned_data = nii_subcube.hdu.data[i, j, k_indices]
aligned_cube = SpectralCube(data=aligned_data,wcs=wcs.WCS(nii_subcube.hdu))




#########################################
##         Main Analysis Part          ##
#########################################

#8um intensity based mask
lower_limits = lower_limits = np.array([ 40.,  58.16861734,  84.58970108, 123.01164881,
       178.8854382 , 260.13796505, 378.2966436 , 550.12481755])
upper_limits = np.array([ 58.16861734,  84.58970108, 123.01164881,
       178.8854382 , 260.13796505, 378.2966436 , 550.12481755,
       800.])
nlimits = len(upper_limits)

m = np.zeros((nlimits,1,1))   # nr. of map copies (1 per threshold)
hdu_dust_copies = m+hdu_dust.data[np.newaxis,:,:]  # create 1 copy of map for each threshold



# Extract velocity values of channels
chan = nii_subcube.spectral_axis
nchan = len(chan)

# make masks
in_feed = (hdu_dust_copies > lower_limits[:,np.newaxis,np.newaxis]) & (hdu_dust_copies < upper_limits[:,np.newaxis,np.newaxis]) # make masks with thresholds

m_in, x_in, y_in = np.where(in_feed) # generate a list of valid positions

#numer of wanted iterations
iterations = 500
n_pixels_in = np.array([np.sqrt(len(m_in[m_in==0])),np.sqrt(len(m_in[m_in==1])),np.sqrt(len(m_in[m_in==2])),\
                           np.sqrt(len(m_in[m_in==3])),np.sqrt(len(m_in[m_in==4])),np.sqrt(len(m_in[m_in==5])),\
                           np.sqrt(len(m_in[m_in==6])),np.sqrt(len(m_in[m_in==7])),np.sqrt(len(m_in[m_in==8]))]).astype(int)

#n_pixels_in = np.array([500,500,500,500,500,500,500,500]).astype(int)

i_in_0 = np.random.randint(0,(len(x_in[m_in<1])),size=(iterations,n_pixels_in[0])) # randomly pick pixel coordinates
i_in_1 = np.random.randint((len(x_in[m_in<1])),(len(x_in[m_in<2])),size=(iterations,n_pixels_in[1])) # randomly pick pixel coordinates
i_in_2 = np.random.randint((len(x_in[m_in<2])),(len(x_in[m_in<3])),size=(iterations,n_pixels_in[2])) # randomly pick pixel coordinates
i_in_3 = np.random.randint((len(x_in[m_in<3])),(len(x_in[m_in<4])),size=(iterations,n_pixels_in[3])) # randomly pick pixel coordinates
i_in_4 = np.random.randint((len(x_in[m_in<4])),(len(x_in[m_in<5])),size=(iterations,n_pixels_in[4])) # randomly pick pixel coordinates
i_in_5 = np.random.randint((len(x_in[m_in<5])),(len(x_in[m_in<6])),size=(iterations,n_pixels_in[5])) # randomly pick pixel coordinates
i_in_6 = np.random.randint((len(x_in[m_in<6])),(len(x_in[m_in<7])),size=(iterations,n_pixels_in[6])) # randomly pick pixel coordinates
i_in_7 = np.random.randint((len(x_in[m_in<7])),(len(x_in[m_in<8])),size=(iterations,n_pixels_in[7])) # randomly pick pixel coordinates


all_spec_in_0 = aligned_cube.hdu.data[:,x_in[i_in_0],y_in[i_in_0]]  # Extract spectra of random pixles
all_spec_in_1 = aligned_cube.hdu.data[:,x_in[i_in_1],y_in[i_in_1]]  # Extract spectra of random pixles
all_spec_in_2 = aligned_cube.hdu.data[:,x_in[i_in_2],y_in[i_in_2]]  # Extract spectra of random pixles
all_spec_in_3 = aligned_cube.hdu.data[:,x_in[i_in_3],y_in[i_in_3]]  # Extract spectra of random pixles
all_spec_in_4 = aligned_cube.hdu.data[:,x_in[i_in_4],y_in[i_in_4]]  # Extract spectra of random pixles
all_spec_in_5 = aligned_cube.hdu.data[:,x_in[i_in_5],y_in[i_in_5]]  # Extract spectra of random pixles
all_spec_in_6 = aligned_cube.hdu.data[:,x_in[i_in_6],y_in[i_in_6]]  # Extract spectra of random pixles
all_spec_in_7 = aligned_cube.hdu.data[:,x_in[i_in_7],y_in[i_in_7]]  # Extract spectra of random pixles

#all_spec_in_norm = all_spec_in/np.max(all_spec_in,axis=0)

av_spec_in = np.zeros((nchan,iterations,8))

av_spec_in[:,:,0] = np.nanmean(all_spec_in_0,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,1] = np.nanmean(all_spec_in_1,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,2] = np.nanmean(all_spec_in_2,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,3] = np.nanmean(all_spec_in_3,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,4] = np.nanmean(all_spec_in_4,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,5] = np.nanmean(all_spec_in_5,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,6] = np.nanmean(all_spec_in_6,axis=2)  # average over the randomly selected pixels
av_spec_in[:,:,7] = np.nanmean(all_spec_in_7,axis=2)  # average over the randomly selected pixels


#
# KS Test on these spectra
#

med_spec_in = np.median(av_spec_in,axis=1) # median spectra

ks_p = np.array([0.,0.,0.,0.,0.,0.,0.])

for i in range(len(ks_p)):
    ks_p[i] = ss.ks_2samp(med_spec_in[:,i],med_spec_in[:,i+1])[1]



#
# Calculating statistics
#       
v = chan.value[:,np.newaxis,np.newaxis]  # velocity
dv = (v[2]-v[1])

P_in = av_spec_in

# Expectation values (E[X])
Ex_in = np.sum(P_in*v*dv,axis=0)/np.sum(P_in*-0.5,axis=0)
Ex_in = Ex_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility

# Std Dev. E[(X-mu)^2]
s_in = np.sqrt(np.sum(((v-Ex_in)**2)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0))  # std
s_in = s_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility

# Skewness E[((X-mu)/sigma)^3]
s3_in = np.sum((((v-Ex_in)/s_in)**3)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0) # moment-3
s3_in = s3_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility

# Kurtosis E[((X-mu)/sigma)^4]
s4_in = np.sum((((v-Ex_in)/s_in)**4)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0) # moment-4
s4_in = s4_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility


std_in = np.median(s_in,axis=(0,1)) # average over extra dim and iterations to obtain value for each threshold
skew_in = np.median(s3_in,axis=(0,1))
kurt_in = np.median(s4_in,axis=(0,1))

err_std_in = np.std(s_in,axis=(0,1))  # standard error in each qnt.
err_skew_in = np.std(s3_in,axis=(0,1))
err_kurt_in = np.std(s4_in,axis=(0,1))

m4_in = np.sum((((v-Ex_in))**4)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0) # moment-4
m4_in = m4_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility
error_std = np.sqrt((((120/121)**2)*(((m4_in/s_in**2)-1)/(4*121))) - (120./(2*(120**3))))
error_in_std = np.mean((error_std)**2,axis=(0,1))

#
# PLotting Section
#

# Plot Spectra

nrows = np.floor(np.sqrt(nlimits)).astype(int)  # Subplot array shape
ncols = np.ceil(nlimits/nrows).astype(int)

width = 5*ncols  # Subplot Size
height = 5*nrows

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(width,height))  # plot initialize
matplotlib.rcParams.update({'font.size': 14})

for ax_x in range(nrows):
    for ax_y in range(ncols):

        
        print((4*ax_x)+ax_y)

        # all spectra (translucent)
        ax[ax_x,ax_y].plot(chan,av_spec_in[:,:,((4*ax_x)+ax_y)],color='b',alpha=0.05)

        ax[ax_x,ax_y].set_title("["+str(lower_limits[(4*ax_x)+ax_y].astype(int))+","+str(upper_limits[(4*ax_x)+ax_y].astype(int))+"] MJy/sr")

        # median spectrum
        ax[ax_x,ax_y].plot(chan,np.median(av_spec_in[:,:,((4*ax_x)+ax_y)],axis=1),color='goldenrod')
        
        ax[ax_x,ax_y].set_ylim(ymax=12)

# Axes Properties
ax[1,0].set_xlabel('Rest frame velocity ($km.s^{-1}$)')
ax[1,0].set_ylabel(r'$T_A^*$ (K)')
#ax[2,0].legend(loc='upper right',fontsize=12)

plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/new_ALP_spectra_intervals.png',dpi=400)



# Statistic vs Threshold Scatter Plot

fig2,ax2 = plt.subplots(nrows=1, ncols=3, figsize=(15,5))  # plot initialize

x_ticks = (upper_limits+lower_limits)/2.0
# Std Dev. Scatter
ax2[0].errorbar(x_ticks,std_in,yerr=err_std_in,color='goldenrod',marker='o',ecolor='blue',mfc='goldenrod',linestyle='none',elinewidth=0.7)
ax2[0].set_xscale('log')
ax2[0].set_xticks(x_ticks)
ax2[0].set_xticklabels(x_ticks.astype(int).astype(str),rotation='vertical')
ax2[0].get_xaxis().set_tick_params(which='minor', size=0)
ax2[0].get_xaxis().set_tick_params(which='minor', width=0) 
ax2[0].set_title("Median Std Dev")  # Axes Properties

# Skewness Scatter
ax2[1].errorbar(x_ticks,skew_in,yerr=err_skew_in,color='goldenrod',marker='o',ecolor='blue',mfc='goldenrod',linestyle='none',elinewidth=0.7)
ax2[1].set_xscale('log')
ax2[1].set_xticks(x_ticks)
ax2[1].set_xticklabels(x_ticks.astype(int).astype(str),rotation='vertical')
ax2[1].get_xaxis().set_tick_params(which='minor', size=0)
ax2[1].get_xaxis().set_tick_params(which='minor', width=0) 
ax2[1].set_title("Median Skewness")

# Kurtosis Scatter
ax2[2].errorbar(x_ticks,kurt_in,yerr=err_kurt_in,color='goldenrod',marker='o',ecolor='blue',mfc='goldenrod',linestyle='none',elinewidth=0.7)
ax2[2].set_xscale('log')
ax2[2].set_xticks(x_ticks)
ax2[2].set_xticklabels(x_ticks.astype(int).astype(str),rotation='vertical')
ax2[2].get_xaxis().set_tick_params(which='minor', size=0)
ax2[2].get_xaxis().set_tick_params(which='minor', width=0) 
ax2[2].set_title("Median Kurtosis")

# Plot Properties
ax2[0].set_xlabel(r'Median 8$\mu m$ Intensity ($MJy.sr^{-1}$)')
ax2[0].set_ylabel('$km.s^{-1}$')
plt.tight_layout()

plt.show()
plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/new_ALP_statistic_intervals.eps')

print(ks_p)
