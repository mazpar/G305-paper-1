#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:30:39 2021

@author: pmazumdar

Task: Multiple analyses based on relation between
      the excitation properties and 8um flux.
"""

#=======================================================================================================================
#                             Plot the 8um Map with 12CO integrated intensity contours                                  
#=======================================================================================================================

import aplpy
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size':12})
fig = plt.figure(figsize=(7,7))
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305-GLM-8-reproject.fits',figure=fig)
F.set_theme('publication')
F.show_colorscale(cmap='Blues',stretch='log',vmin=100)
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_font(size='medium', weight='medium', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')
F.colorbar.set_location('right')
F.colorbar.set_width(0.1)
F.colorbar.set_pad(0.07)
F.colorbar.set_axis_label_text(r"(MJy/Sr)")
#F.colorbar.set_axis_label_pad(30)
#F.colorbar.set_axis_label_rotation(0)
#F.axis_labels.hide()
#F.set_title(r"Excitation Temperature (K)")

#line_list = [np.array([[305.351,305.0788],[0.1210355,0.4414058]]),np.array([[305.3177,305.2177],[0.009924373,-0.1789645]])]
#F.show_lines(line_list,color='white',linewidth=1.2,linestyle='dashed')
#F.add_label(305.2,0.33,'N',color='white',family='serif',weight='bold',size=12)
#F.add_label(305.28,-0.1,'S',color='white',family='serif',weight='bold',size=12)
#F.show_lines(line_list,color='white',linewidth=1.2,linestyle='dashed')
#F.add_label(305.2,0.33,'N',color='white',family='serif',weight='bold',size=14)
#F.add_label(305.28,0.06,'S',color='white',family='serif',weight='bold',size=14)

#level=np.array([30,70,120])
level=np.array([4,20,80])
#F.show_contour(data='/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO-moment0.fits', \
#               levels=level, colors='k', smooth=1, linewidths=0.8, returnlevels=True)
F.show_contour(data='/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO-moment0.fits', \
               levels=level, smooth=1, colors='k', linewidths=0.8, returnlevels=True)

F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/glimpse_8u_overlay_13CO.eps')


#=======================================================================================================================
#                       Histogram Distribution of Integrated properties vs 8um Flux
#=======================================================================================================================


from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
import matplotlib
import aplpy
from scipy.stats import median_absolute_deviation as mad
from scipy.stats import gaussian_kde

# Read Ratio File
hdu_ratio = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305-excitation.fits')[0]
hdu_temp =fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/tex_12_mean_squeeze.fits')[0]
hdu_cd =fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal_sum_new.fits')[0]

# Compare excitation with 8um
hdu_dust = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
w_dust = wcs.WCS(hdu_dust.header)


flat_dust_ratio = hdu_dust.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_ratio.data.flatten()))]
flat_ratio = hdu_ratio.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_ratio.data.flatten()))]
flat_dust_temp = hdu_dust.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_temp.data.flatten()))& (~np.isnan(hdu_ratio.data.flatten()))]
flat_temp = hdu_temp.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_temp.data.flatten()))& (~np.isnan(hdu_ratio.data.flatten()))]
flat_dust_temp_nan = hdu_dust.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_temp.data.flatten()))& (np.isnan(hdu_ratio.data.flatten()))]
flat_temp_nan = hdu_temp.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_temp.data.flatten()))& (np.isnan(hdu_ratio.data.flatten()))]
flat_dust_cd = hdu_dust.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_cd.data.flatten()))]
flat_cd = hdu_cd.data.flatten()[~np.isnan(hdu_dust.data.flatten())& (~np.isnan(hdu_cd.data.flatten()))]


# =============================================================================
# Plot Histograms
# =============================================================================

limits = np.logspace(np.log10(40),np.log10(800),9)

#plt.rcParams.update({'font.size': 14})
#from astropy.visualization import hist

# =============================================================================
# T_ex (Ratio) Histogram
# =============================================================================

#fig,ax = plt.subplots(1, 2, figsize=(12, 6))
#fig.subplots_adjust(bottom=0.15)

#limits = [60,80,100,125,150,200,250,300,400]
meds = [0,0,0,0,0,0,0,0]
q25 = [0,0,0,0,0,0,0,0]
q75 = [0,0,0,0,0,0,0,0]
err = [0,0,0,0,0,0,0,0]
med_8um_flux = [0,0,0,0,0,0,0,0]
for i in range(len(limits)-1):
    flat_ratio_interval = flat_ratio[(flat_dust_ratio>limits[i]) & (flat_dust_ratio<limits[i+1])]
    meds[i] = np.median(flat_ratio_interval)
    q25[i] = np.percentile(flat_ratio_interval,25)
    q75[i] = np.percentile(flat_ratio_interval,75)    
#    err[i] = mad(flat_temp_interval)
    med_8um_flux[i] = 10**((np.log10(limits[i])+np.log10(limits[i+1]))/2.0)
#    hist(flat_ratio_interval,bins='knuth',histtype='step', ax=ax[0], alpha=0.4, density=1, color=plt.cm.Reds((1-i/(1.2*len(limits)))), \
#         label=str(limits[i].astype(int))+"-"+str(limits[i+1].astype(int)))

#ax[0].set_xlabel(r'Ratio $^{13}$CO (3-2)/(2-1)')
#ax[0].set_ylabel('Pixel Density')
#ax[0].set_xlim(xmax=1.25)
#ax[0].legend()

#ax[1].errorbar(med_8um_flux,meds,yerr=err,color='r',ecolor='k',marker='s',linestyle='none')
#ax[1].scatter(med_8um_flux,meds,color='r',marker='s')
#ax[1].vlines(med_8um_flux, q25, q75, color='k', linestyle='-', lw=1)
#ax[1].set_xscale('log')
#ax[1].set_xlim(xmax=1000)
#ax[1].set_xlabel(r'Median of 8$\mu$m Interval')
#ax[1].set_ylabel(r'Median $^{13}$CO (3-2)/(2-1) Ratio')
#plt.tight_layout()

#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/excitation-vs-8um.eps",dvi=300)

# =============================================================================
# Ratio Scatter
# =============================================================================
xy = np.vstack([flat_dust_ratio,flat_ratio])
z = gaussian_kde(xy)(xy)
matplotlib.rcParams.update({'font.size': 14})
fig,ax=plt.subplots(1,1,figsize=(7,6))
plt.scatter(flat_dust_ratio,flat_ratio,s=5,c=z,cmap=plt.cm.nipy_spectral)
cbar = plt.colorbar(aspect=30,pad=0.01)
colors = ['k','k','k','k','k','w','w','m']
plt.scatter(med_8um_flux,meds,color=colors,marker='o')
plt.vlines(med_8um_flux, q25, q75, color=colors, linestyle='-', lw=2)
plt.xscale('log')
plt.ylim(top=1.6)
plt.xlabel('8$\mu$m Flux [MJy.sr$^{-1}$]')
plt.ylabel('$^{13}$CO (3-2)/(2-1) Ratio')
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/ratio-vs-8um-scatter.eps",dvi=60)

# =============================================================================
# Temp Histogram
# =============================================================================

#fig,ax = plt.subplots(1, 2, figsize=(12, 6))
#fig.subplots_adjust(bottom=0.15)

#limits = [60,80,100,125,150,200,250,300,400]
meds = [0,0,0,0,0,0,0,0]
err = [0,0,0,0,0,0,0,0]
q25 = [0,0,0,0,0,0,0,0]
q75 = [0,0,0,0,0,0,0,0]
med_8um_flux = [0,0,0,0,0,0,0,0]
for i in range(len(limits)-1):
  flat_temp_interval = flat_temp[(flat_dust_temp>limits[i]) & (flat_dust_temp<limits[i+1])]
  meds[i] = np.median(flat_temp_interval)
  q25[i] = np.percentile(flat_temp_interval,25)
  q75[i] = np.percentile(flat_temp_interval,75)    
#    err[i] = mad(flat_temp_interval)
  med_8um_flux[i] = 10**((np.log10(limits[i])+np.log10(limits[i+1]))/2.0)
#    hist(flat_temp_interval,bins='knuth',histtype='step', ax=ax[0], alpha=0.4, density=1, color=plt.cm.Reds((1-i/(1.2*len(limits)))), \
#         label=str(limits[i].astype(int))+"-"+str(limits[i+1].astype(int)))

#ax[0].set_xlabel(r'Excitation Temperature (K)')
#ax[0].set_ylabel('Pixel Density')
#ax[0].legend()

#ax[1].errorbar(med_8um_flux,meds,yerr=err,color='r',ecolor='k',marker='s',linestyle='none')
#ax[1].scatter(med_8um_flux,meds,color='r',marker='s')
#ax[1].vlines(med_8um_flux, q25, q75, color='k', linestyle='-', lw=1)
#ax[1].set_xscale('log')
#ax[1].set_xlim(xmax=1000)
#ax[1].set_xlabel(r'Median of 8$\mu$m Interval')
#ax[1].set_ylabel(r'Median Excitation Temperature (K)')
#plt.tight_layout()

#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Temperature-vs-8um.eps",dvi=300)

# =============================================================================
# Temp Scatter
# =============================================================================
matplotlib.rcParams.update({'font.size': 14})
xy = np.vstack([flat_dust_temp,flat_temp])
z = gaussian_kde(xy)(xy)
fig,ax=plt.subplots(1,1,figsize=(7,6))
plt.scatter(flat_dust_temp,flat_temp,s=5,c=z,cmap=plt.cm.nipy_spectral)
cbar = plt.colorbar(aspect=30,pad=0.01)
plt.scatter(flat_dust_temp_nan,flat_temp_nan,s=2,c='gray',alpha=0.05)
colors = ['k','k','w','w','w','w','w','m']
plt.scatter(med_8um_flux,meds,color=colors,marker='o')
plt.vlines(med_8um_flux, q25, q75, color=colors, linestyle='-', lw=2)

def pow_func(x,a,b):
  return a*(x**b)

import scipy.optimize as so
popt, pcov = so.curve_fit(pow_func, flat_dust_temp, flat_temp)
perr = np.sqrt(np.diag(pcov))
xdata = np.logspace(1, 4, 50)
plt.plot(xdata, pow_func(xdata, *popt), 'b',label=r'Fit: %5.3f $\pm$ %5.3f x$^{%5.3f \pm %5.3f}$' %(popt[0],perr[0],popt[1],perr[1]))
plt.legend(loc='upper right')

plt.xscale('log')
plt.yscale('log')
plt.xlim(left=20)
plt.yticks([10,20,30,40])

plt.xlabel('8$\mu$m Flux [MJy.sr$^{-1}$]')
plt.ylabel('Excitation Temperature [K]')

plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/excitation-vs-8um-scatter.eps",dvi=60)



# =============================================================================
# Column Density Histogram
# =============================================================================

#fig,ax = plt.subplots(1, 2, figsize=(12, 6))
#fig.subplots_adjust(bottom=0.15)

#limits = [60,80,100,125,150,200,250,300,400]
meds = [0,0,0,0,0,0,0,0]
#err = [0,0,0,0,0,0,0,0]
q25 = [0,0,0,0,0,0,0,0]
q75 = [0,0,0,0,0,0,0,0]
med_8um_flux = [0,0,0,0,0,0,0,0]
cd_interval = [[],[],[],[],[],[],[],[]]
for i in range(len(limits)-1):
  cd_interval[i] = flat_cd[(flat_dust_cd>limits[i]) & (flat_dust_cd<limits[i+1])]
  flat_cd_interval = flat_cd[(flat_dust_cd>limits[i]) & (flat_dust_cd<limits[i+1])]
  meds[i] = np.median(flat_cd_interval)
  q25[i] = np.percentile(flat_cd_interval,25)
  q75[i] = np.percentile(flat_cd_interval,75)    
#	err[i] = mad(flat_cd_interval)
  med_8um_flux[i] = 10**((np.log10(limits[i])+np.log10(limits[i+1]))/2.0)
#	hist(flat_cd_interval,bins='knuth',histtype='step', ax=ax[0], alpha=0.4, density=1, color=plt.cm.Reds((1-i/(1.2*len(limits)))), \
#         label=str(limits[i].astype(int))+"-"+str(limits[i+1].astype(int)))

#ax[0].set_ylim(ymax=0.9e-16)
#ax[0].set_xscale('log')
#ax[0].set_xlabel(r'$^{13}$CO Column Density (cm$^{-2}$)')
#ax[0].set_ylabel('Pixel Density')
#ax[0].legend()



# Statistic vs Threshold Violin Plot


#def set_axis_style(ax, labels):  # function for axes properties
#    ax.get_xaxis().set_tick_params(direction='out')
#    ax.xaxis.set_ticks_position('bottom')
#    ax.set_xticks(np.arange(1,len(labels)+1))
#    ax.set_xticklabels(np.floor(labels).astype(str), rotation=90, ha='center')
    

#labels = np.array(med_8um_flux) # Custom x-Labels for plot

#violins = ax[1].violinplot(cd_interval,showmedians=False,showextrema=False)

#quartile1 = []
#median = []
#quartile3 = []

#for i in range(len(cd_interval)):
#	quartile1 = np.append(quartile1,np.percentile(cd_interval[i], 25))
#	quartile3 = np.append(quartile3,np.percentile(cd_interval[i], 75))
#	median = np.append(median,np.percentile(cd_interval[i], 50))
	 	 
#inds = np.arange(1, len(labels) + 1)
#ax[1].scatter(inds, median, marker='o', color='k', s=30)
#ax[1].vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=1)


#for violin in violins['bodies']:
#        violin.set_facecolor('mediumseagreen')
#        violin.set_edgecolor('black')
#        violin.set_alpha(0.4)

# general axes properties
#set_axis_style(ax[1], labels)
#ax[1].set_yscale('log')
#ax[1].set_xlabel(r'8$\mu m$ median flux ($MJy.sr^{-1}$)')
#ax[1].set_ylabel('Column Density (cm$^{-2}$)')

#plt.tight_layout()
#plt.show()

#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Column-density-vs-8um.eps",dvi=300)

# =============================================================================
# Column Density Scatter
# =============================================================================
matplotlib.rcParams.update({'font.size': 12})
xy = np.vstack([flat_dust_cd,flat_cd])
z = gaussian_kde(xy)(xy)
fig,ax=plt.subplots(1,1,figsize=(7,6))
plt.scatter(flat_dust_cd,flat_cd,s=5,c=z,cmap=plt.cm.nipy_spectral)
cbar = plt.colorbar(aspect=30,pad=0.01)
colors = ['k','k','k','k','k','w','w','m']
plt.scatter(med_8um_flux,meds,color=colors,marker='o')
plt.vlines(med_8um_flux, q25, q75, color=colors, linestyle='-', lw=2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('8$\mu$m Flux [MJy.sr$^{-1}$]')
plt.ylabel('Column Density (cm$^{-2}$)')
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/ncd-vs-8um-scatter.eps",dvi=60)




#=======================================================================================================================
#                      8um Map with contours of levels
#=======================================================================================================================
import matplotlib.markers as markers

fig = plt.figure(figsize=(12,12))
matplotlib.rcParams.update({'font.size': 9})
r = plt.cm.Reds
#levs = limits
#levs = [84.0,123.0,178.0,260.,500.]
levs = [60,150]
colors = ['tab:purple','tab:brown','tab:blue','tab:green','tab:orange']
#F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305-GLM-8-reproject.fits',figure=fig)
#F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO-moment0.fits',figure=fig)
#F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/tex_12_mean_squeeze.fits',figure=fig)
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal_sum_new.fits',figure=fig)
F.set_theme('publication')
F.show_grayscale(stretch='log')
F.set_nan_color('white')
#F.add_colorbar()
le = F.show_contour(data = hdu_dust.data,levels=levs,colors=colors,linewidths=0.8,smooth=3,returnlevels=1)

xrms = np.array([305.6812])
yrms = np.array([0.0985])
rms_marker = markers.MarkerStyle(marker=(4,0,90),fillstyle="none")
F.show_markers(xrms,yrms,edgecolor='red',facecolor='none',marker=rms_marker,s=20)

xuchii = np.array([305.362,305.368,305.562,305.55,305.200])
yuchii = np.array([0.150,0.213,0.013,-0.012,0.019])
uchii_marker = markers.MarkerStyle(marker=(4,0,0),fillstyle="none")
F.show_markers(xuchii,yuchii,edgecolor='red',facecolor='none',marker=uchii_marker,s=20)

xhii = np.array([305.353,305.195,305.270,305.320,305.254,305.348,305.551,305.532])
yhii = np.array([0.193,0.033,-0.007,0.070,0.204,0.223,0.014,0.348])
hii_marker = markers.MarkerStyle(marker=(3,0,0),fillstyle="none")
F.show_markers(xhii,yhii,facecolor='none',edgecolor='cyan',marker=hii_marker,s=20)

F.save("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_13CO_8um_feedback_contour.eps")
#le = F.show_contour(data = hdu_dust.data , cmap=plt.cm.Reds, levels=levs,linewidths=0.8,smooth=1,returnlevels=1)
#F.add_label(0.5, 0.95, 'levels='+str(limits)+ 'MJy/sr', color='white',relative=True)
#F.save("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/8um-contours.eps")


