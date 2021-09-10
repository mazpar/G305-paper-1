#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:17:58 2020

@author: pmazumdar

Task: read Paper-1 section 7.2

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

#from scipy.optimize import curve_fit
#from matplotlib.ticker import ScalarFormatter
#from astropy.coordinates import Angle, SkyCoord
#from regions import CircleSkyRegion
#from regions.core import PixCoord
#from regions.shapes.circle import CirclePixelRegion


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
#mask_12CO = mask_nan & mask_uchii
mask_12CO = mask_nan

nii_subcube = nii_subcube.with_mask(mask_12CO)   # apply mask to spectral cube




#################################################################
# Aligning spectra based on moment 1 map / from velocity peaks ##
#################################################################

#----> mom-1 based shifts <----
mom1_data = hdu_mom1.data
# mom1_data[np.isnan(mom1_data)] = -40
mom1_data[mom1_data<-50] = -40 # values <-60 are pulled back to -40 since they are most likely from noise
mom1_data[mom1_data>-15] = -40 # values >-10 are pulled back to -40 since they are most likely from noise

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
#mask_limits = np.array([40,80,120,160,200,250,350,450,600])
mask_limits = np.logspace(np.log10(40),np.log10(800),9)
nlimits = len(mask_limits)
m = np.zeros((nlimits,1,1))   # nr. of map copies (1 per threshold)
hdu_dust_copies = m+hdu_dust.data[np.newaxis,:,:]  # create 1 copy of map for each threshold


# Extract velocity values of channels
chan = nii_subcube.spectral_axis
nchan = len(chan)

# make masks
in_feed = hdu_dust_copies > mask_limits[:,np.newaxis,np.newaxis] # make masks with thresholds

m_in, x_in, y_in = np.where(in_feed) # generate a list of valid positions
m_out, x_out, y_out = np.where(~in_feed) # generate a list of valid positions

#numer of wanted iterations
iterations = 1000
n_pixels_in = np.array([np.sqrt(len(m_in[m_in==0])),np.sqrt(len(m_in[m_in==1])),np.sqrt(len(m_in[m_in==2])),\
                           np.sqrt(len(m_in[m_in==3])),np.sqrt(len(m_in[m_in==4])),np.sqrt(len(m_in[m_in==5])),\
                           np.sqrt(len(m_in[m_in==6])),np.sqrt(len(m_in[m_in==7])),np.sqrt(len(m_in[m_in==8]))]).astype(int)
n_pixels_out = np.array([np.sqrt(len(m_out[m_out==0])),np.sqrt(len(m_out[m_out==1])),np.sqrt(len(m_out[m_out==2])),\
                           np.sqrt(len(m_out[m_out==3])),np.sqrt(len(m_out[m_out==4])),np.sqrt(len(m_out[m_out==5])),\
                           np.sqrt(len(m_out[m_out==6])),np.sqrt(len(m_out[m_out==7])),np.sqrt(len(m_out[m_out==8]))]).astype(int)

i_in_0 = np.random.randint(0,(len(x_in[m_in<1])),size=(iterations,n_pixels_in[0])) # randomly pick pixel coordinates
i_out_0 = np.random.randint(0,(len(x_out[m_out<1])),size=(iterations,n_pixels_out[0])) # randomly pick pixel coordinates
i_in_1 = np.random.randint((len(x_in[m_in<1])),(len(x_in[m_in<2])),size=(iterations,n_pixels_in[1])) # randomly pick pixel coordinates
i_out_1 = np.random.randint((len(x_out[m_out<1])),(len(x_out[m_out<2])),size=(iterations,n_pixels_out[1])) # randomly pick pixel coordinates
i_in_2 = np.random.randint((len(x_in[m_in<2])),(len(x_in[m_in<3])),size=(iterations,n_pixels_in[2])) # randomly pick pixel coordinates
i_out_2 = np.random.randint((len(x_out[m_out<2])),(len(x_out[m_out<3])),size=(iterations,n_pixels_out[2])) # randomly pick pixel coordinates
i_in_3 = np.random.randint((len(x_in[m_in<3])),(len(x_in[m_in<4])),size=(iterations,n_pixels_in[3])) # randomly pick pixel coordinates
i_out_3 = np.random.randint((len(x_out[m_out<3])),(len(x_out[m_out<4])),size=(iterations,n_pixels_out[3])) # randomly pick pixel coordinates
i_in_4 = np.random.randint((len(x_in[m_in<4])),(len(x_in[m_in<5])),size=(iterations,n_pixels_in[4])) # randomly pick pixel coordinates
i_out_4 = np.random.randint((len(x_out[m_out<4])),(len(x_out[m_out<5])),size=(iterations,n_pixels_out[4])) # randomly pick pixel coordinates
i_in_5 = np.random.randint((len(x_in[m_in<5])),(len(x_in[m_in<6])),size=(iterations,n_pixels_in[5])) # randomly pick pixel coordinates
i_out_5 = np.random.randint((len(x_out[m_out<5])),(len(x_out[m_out<6])),size=(iterations,n_pixels_out[5])) # randomly pick pixel coordinates
i_in_6 = np.random.randint((len(x_in[m_in<6])),(len(x_in[m_in<7])),size=(iterations,n_pixels_in[6])) # randomly pick pixel coordinates
i_out_6 = np.random.randint((len(x_out[m_out<6])),(len(x_out[m_out<7])),size=(iterations,n_pixels_out[6])) # randomly pick pixel coordinates
i_in_7 = np.random.randint((len(x_in[m_in<7])),(len(x_in[m_in<8])),size=(iterations,n_pixels_in[7])) # randomly pick pixel coordinates
i_out_7 = np.random.randint((len(x_out[m_out<7])),(len(x_out[m_out<8])),size=(iterations,n_pixels_out[7])) # randomly pick pixel coordinates
i_in_8 = np.random.randint((len(x_in[m_in<8])),(len(x_in[m_in<9])),size=(iterations,n_pixels_in[8])) # randomly pick pixel coordinates
i_out_8 = np.random.randint((len(x_out[m_out<8])),(len(x_out[m_out<9])),size=(iterations,n_pixels_out[8])) # randomly pick pixel coordinates


all_spec_in_0 = aligned_cube.hdu.data[:,x_in[i_in_0],y_in[i_in_0]]  # Extract spectra of random pixles
all_spec_out_0 = aligned_cube.hdu.data[:,x_out[i_out_0],y_out[i_out_0]]
all_spec_in_1 = aligned_cube.hdu.data[:,x_in[i_in_1],y_in[i_in_1]]  # Extract spectra of random pixles
all_spec_out_1 = aligned_cube.hdu.data[:,x_out[i_out_1],y_out[i_out_1]]
all_spec_in_2 = aligned_cube.hdu.data[:,x_in[i_in_2],y_in[i_in_2]]  # Extract spectra of random pixles
all_spec_out_2 = aligned_cube.hdu.data[:,x_out[i_out_2],y_out[i_out_2]]
all_spec_in_3 = aligned_cube.hdu.data[:,x_in[i_in_3],y_in[i_in_3]]  # Extract spectra of random pixles
all_spec_out_3 = aligned_cube.hdu.data[:,x_out[i_out_3],y_out[i_out_3]]
all_spec_in_4 = aligned_cube.hdu.data[:,x_in[i_in_4],y_in[i_in_4]]  # Extract spectra of random pixles
all_spec_out_4 = aligned_cube.hdu.data[:,x_out[i_out_4],y_out[i_out_4]]
all_spec_in_5 = aligned_cube.hdu.data[:,x_in[i_in_5],y_in[i_in_5]]  # Extract spectra of random pixles
all_spec_out_5 = aligned_cube.hdu.data[:,x_out[i_out_5],y_out[i_out_5]]
all_spec_in_6 = aligned_cube.hdu.data[:,x_in[i_in_6],y_in[i_in_6]]  # Extract spectra of random pixles
all_spec_out_6 = aligned_cube.hdu.data[:,x_out[i_out_6],y_out[i_out_6]]
all_spec_in_7 = aligned_cube.hdu.data[:,x_in[i_in_7],y_in[i_in_7]]  # Extract spectra of random pixles
all_spec_out_7 = aligned_cube.hdu.data[:,x_out[i_out_7],y_out[i_out_7]]
all_spec_in_8 = aligned_cube.hdu.data[:,x_in[i_in_8],y_in[i_in_8]]  # Extract spectra of random pixles
all_spec_out_8 = aligned_cube.hdu.data[:,x_out[i_out_8],y_out[i_out_8]]

#all_spec_in[0:40,:,:,:] = 0.0
#all_spec_out[0:40,:,:,:] = 0.0

#all_spec_in_norm = all_spec_in/np.max(all_spec_in,axis=0)
#all_spec_out_norm = all_spec_out/np.max(all_spec_out,axis=0)

av_spec_in = np.zeros((nchan,iterations,9))
av_spec_out = np.zeros((nchan,iterations,9))

av_spec_in[:,:,0] = np.nanmean(all_spec_in_0,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,0] = np.nanmean(all_spec_out_0,axis=2)
av_spec_in[:,:,1] = np.nanmean(all_spec_in_1,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,1] = np.nanmean(all_spec_out_1,axis=2)
av_spec_in[:,:,2] = np.nanmean(all_spec_in_2,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,2] = np.nanmean(all_spec_out_2,axis=2)
av_spec_in[:,:,3] = np.nanmean(all_spec_in_3,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,3] = np.nanmean(all_spec_out_3,axis=2)
av_spec_in[:,:,4] = np.nanmean(all_spec_in_4,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,4] = np.nanmean(all_spec_out_4,axis=2)
av_spec_in[:,:,5] = np.nanmean(all_spec_in_5,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,5] = np.nanmean(all_spec_out_5,axis=2)
av_spec_in[:,:,6] = np.nanmean(all_spec_in_6,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,6] = np.nanmean(all_spec_out_6,axis=2)
av_spec_in[:,:,7] = np.nanmean(all_spec_in_7,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,7] = np.nanmean(all_spec_out_7,axis=2)
av_spec_in[:,:,8] = np.nanmean(all_spec_in_8,axis=2)  # average over the randomly selected pixels
av_spec_out[:,:,8] = np.nanmean(all_spec_out_8,axis=2)




#av_spec_in = np.nanmean(all_spec_in_norm,axis=2)  # average over the randomly selected pixels
#av_spec_out = np.nanmean(all_spec_out_norm,axis=2)

med_spec_in = np.median(av_spec_in,axis=1) # median spectra
med_spec_out = np.median(av_spec_out,axis=1) # median spectra

ks_in = np.array([0.,0.,0.,0.,0.,0.,0.,0.])
ks_out = np.array([0.,0.,0.,0.,0.,0.,0.,0.])
ks_cross = np.array([0.,0.,0.,0.,0.,0.,0.,0.])

for i in range(len(ks_in)):
    ks_in[i] = ss.ks_2samp(med_spec_in[:,i],med_spec_in[:,i+1])[1]
    ks_out[i] = ss.ks_2samp(med_spec_out[:,i],med_spec_out[:,i+1])[1]
    ks_cross[i] = ss.ks_2samp(med_spec_out[:,i],med_spec_in[:,i])[1]
#
# Calculating statistics
#       
v = chan.value[:,np.newaxis,np.newaxis]  # velocity
dv = (v[2]-v[1])

P_in = av_spec_in
P_out = av_spec_out

# Expectation values (E[X])
Ex_in = np.sum(P_in*v*dv,axis=0)/np.sum(P_in*-0.5,axis=0)
Ex_out =  np.sum(P_out*v*dv,axis=0)/np.sum(P_out*dv,axis=0)
Ex_in = Ex_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility
Ex_out = Ex_out[np.newaxis,:,:]

# Std Dev. E[(X-mu)^2]
s_in = np.sqrt(np.sum(((v-Ex_in)**2)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0))  # std
s_out = np.sqrt(np.sum(((v-Ex_out)**2)*P_out*dv,axis=0)/np.sum(P_out*dv,axis=0))
s_in = s_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility
s_out = s_out[np.newaxis,:,:]

# Skewness E[((X-mu)/sigma)^3]
s3_in = np.sum((((v-Ex_in)/s_in)**3)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0) # moment-3
s3_out = np.sum((((v-Ex_out)/s_out)**3)*P_out*dv,axis=0)/np.sum(P_out*dv,axis=0)
s3_in = s3_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility
s3_out = s3_out[np.newaxis,:,:]

# Kurtosis E[((X-mu)/sigma)^4]
s4_in = np.sum((((v-Ex_in)/s_in)**4)*P_in*dv,axis=0)/np.sum(P_in*dv,axis=0) # moment-4
s4_out = np.sum((((v-Ex_out)/s_out)**4)*P_out*dv,axis=0)/np.sum(P_out*dv,axis=0)
s4_in = s4_in[np.newaxis,:,:]  # add dimension for arithmetic compatibility
s4_out = s4_out[np.newaxis,:,:]


std_in = np.median(s_in,axis=(0,1)) # average over extra dim and iterations to obtain value for each threshold
std_out = np.median(s_out,axis=(0,1))
skew_in = np.median(s3_in,axis=(0,1))
skew_out = np.median(s3_out,axis=(0,1))
kurt_in = np.median(s4_in,axis=(0,1))
kurt_out = np.median(s4_out,axis=(0,1))

err_std_in = np.std(s_in,axis=(0,1))  # standard error in each qnt.
err_std_out = np.std(s_out,axis=(0,1))
err_skew_in = np.std(s3_in,axis=(0,1))
err_skew_out = np.std(s3_out,axis=(0,1))
err_kurt_in = np.std(s4_in,axis=(0,1))
err_kurt_out = np.std(s4_out,axis=(0,1))

# Finding channel number of peak T
ref_med_in = np.median(av_spec_in[:,:,0],axis=1)
ref_med_out = np.median(av_spec_out[:,:,0],axis=1)
#
# PLotting Section
#

# Plot Spectra

ncols = np.floor(np.sqrt(nlimits)).astype(int)  # Subplot array shape
nrows = np.ceil(nlimits/ncols).astype(int)

width = 5*ncols  # Subplot Size
height = 5*nrows

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(width,height))  # plot initialize
matplotlib.rcParams.update({'font.size': 14})

for ax_x in range(nrows):
    for ax_y in range(ncols):

        print((3*ax_x)+ax_y)

        # all spectra (translucent)
        ax[ax_x,ax_y].plot(chan,av_spec_in[:,:,((3*ax_x)+ax_y)],color='r',alpha=0.05)
        ax[ax_x,ax_y].plot(chan,av_spec_out[:,:,((3*ax_x)+ax_y)],color='b',alpha=0.05)

        # median spectrum
        ax[ax_x,ax_y].plot(chan,np.median(av_spec_in[:,:,((3*ax_x)+ax_y)],axis=1),color='k',label="Inside Feedback Zone")
        ax[ax_x,ax_y].plot(chan,np.median(av_spec_out[:,:,((3*ax_x)+ax_y)],axis=1),color='goldenrod',label="Outside Feedback Zone")
        ax[ax_x,ax_y].axvline(-0.5*np.where(np.median(av_spec_in[:,:,0],axis=1)==(np.max(np.median(av_spec_in[:,:,0],axis=1))))[0],color='k',linestyle='dotted')
        ax[ax_x,ax_y].axvline(-0.5*np.where(np.median(av_spec_out[:,:,0],axis=1)==(np.max(np.median(av_spec_out[:,:,0],axis=1))))[0],color='goldenrod',linestyle='dotted')
        ax[ax_x,ax_y].set_title(str(np.floor(mask_limits[(3*ax_x)+ax_y]))+" MJy/sr")

# Axes Properties
ax[2,0].set_xlabel('Rest frame velocity ($km.s^{-1}$)')
ax[2,0].set_ylabel(r'$T_A^*$ (K)')
ax[2,0].legend(loc='upper right',fontsize=11,fancybox=True, framealpha=0.5)

plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/new_ALP_spectra_threshold.png',dpi=400)



# Statistic vs Threshold Violin Plot

fig3,ax3 = plt.subplots(nrows=1, ncols=3, figsize=(15,5)) # plot initialize


def set_axis_style(ax, labels, title):  # function for axes properties
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=90, ha='center')
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_title(title)
    

labels = mask_limits.astype(int).astype(str)  # Custom x-Labels for plot
titles = ["Standard Deviation","Skewness","Kurtosis"]

# Load Data
data_in = [s_in.reshape(1000,9),s3_in.reshape(1000,9),s4_in.reshape(1000,9)]
data_out = [s_out.reshape(1000,9),s3_out.reshape(1000,9),s4_out.reshape(1000,9)]

# plotting
for a, axis in enumerate(ax3) :
    #parts_in = axis.violinplot(data_in[a],showmedians=False,showextrema=False)
    #parts_out = axis.violinplot(data_out[a],showmedians=False,showextrema=False)

    # appearence of violin plots
    #for pc_in in parts_in['bodies']:
    #    pc_in.set_facecolor('red')
    #    pc_in.set_edgecolor('red')
    #    pc_in.set_alpha(0.0)

    #for pc_out in parts_out['bodies']:
    #    pc_out.set_facecolor('blue')
    #    pc_out.set_edgecolor('blue')
    #    pc_out.set_alpha(0.0)
    
    # Quartile Values of data to plot
    quartile1_in, medians_in, quartile3_in = np.percentile(data_in[a], [25, 50, 75], axis=0)
    quartile1_out, medians_out, quartile3_out = np.percentile(data_out[a], [25, 50, 75], axis=0)

    # plot quartiles as bars over violin plots
    inds = np.arange(1, len(medians_in) + 1)
    axis.scatter(inds, medians_in, marker='o', color='k', s=30, zorder=3,label="Inside Feedback Zone")
    axis.vlines(inds, quartile1_in, quartile3_in, color='red', linestyle='-', lw=0.7)
    axis.scatter(inds, medians_out, marker='o', color='goldenrod', s=30, zorder=3,label="Outside Feedback Zone")
    axis.vlines(inds, quartile1_out, quartile3_out, color='blue', linestyle='-', lw=0.7)
    
    # set axis style
    set_axis_style(axis, labels, titles[a])

# general axes properties
ax3[0].set_xlabel(r'8$\mu m$ threshold ($MJy.sr^{-1}$)')
ax3[0].set_ylabel('$km.s^{-1}$')
ax3[0].legend()
plt.tight_layout()
plt.show()
plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/new_ALP_statistic_threshold.eps')
