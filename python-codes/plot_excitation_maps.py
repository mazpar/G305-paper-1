"""
Created on Tue Aug 25 15:03:59 2020

@author: pmazumdar

Task: Plotting Excitation Maps
"""

import aplpy
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size':16})
fig = plt.figure(figsize=(12,12))
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/tex_12_mean_squeeze.fits',figure=fig)
F.set_theme('publication')
F.show_colorscale(cmap='coolwarm')
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_font(size='medium', weight='medium', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')
F.colorbar.set_location('right')
F.colorbar.set_width(0.14)
F.colorbar.set_pad(0.07)
F.colorbar.set_axis_label_text(r"[K]")
F.colorbar.set_axis_label_pad(30)
F.colorbar.set_axis_label_rotation(0)
#F.axis_labels.hide()
F.set_title(r"Excitation Temperature (K)")

import matplotlib.markers as markers
xhii = np.array([305.353,305.195,305.270,305.320,305.254,305.348,305.551,305.532])
yhii = np.array([0.193,0.033,-0.007,0.070,0.204,0.223,0.014,0.348])
hii_marker = markers.MarkerStyle(marker=(3,0,0),fillstyle="none")
F.show_markers(xhii,yhii,edgecolor='k',facecolor='none',marker=hii_marker,s=70)

xuchii = np.array([305.362,305.368,305.562,305.55,305.200])
yuchii = np.array([0.150,0.213,0.013,-0.012,0.019])
uchii_marker = markers.MarkerStyle(marker=(4,0,0),fillstyle="none")
F.show_markers(xuchii,yuchii,c='k',marker=uchii_marker,s=70)

xbrc = np.array([305.244])
ybrc = np.array([0.224])
brc_marker = markers.MarkerStyle(marker=(4,1,0),fillstyle="none")
F.show_markers(xbrc,ybrc,marker=brc_marker,s=70,mfc='none')

line_list = [np.array([[305.351,305.0788],[0.1210355,0.4414058]]),np.array([[305.3177,305.2177],[0.009924373,-0.1789645]])]
F.show_lines(line_list,color=plt.cm.gray(0.3),linewidth=1.2,linestyle='dashed')
F.add_label(305.2,0.33,'N',color=plt.cm.gray(0.3),family='serif',weight='bold',size=14)
F.add_label(305.28,-0.1,'S',color=plt.cm.gray(0.3),family='serif',weight='bold',size=14)
F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/tex_12_mean_BIG.eps')

plt.rcParams.update({'font.size':16})
fig1 = plt.figure(figsize=(12,12))
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal_sum_new.fits',figure=fig1)
F.set_theme('publication')
F.show_colorscale(cmap='cividis',stretch='log')
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_font(size='medium', weight='medium', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')
F.colorbar.set_location('right')
F.colorbar.set_width(0.14)
F.colorbar.set_pad(0.07)
#F.colorbar.set_axis_label_text(r"$^{13}$CO Column Density (cm$^{-2}$)")
#F.colorbar.set_axis_label_pad(30)
#F.colorbar.set_axis_label_rotation(0)
F.set_title(r"$^{13}$CO Column Density (cm$^{-2}$)")
#F.axis_labels.hide()
F.show_lines(line_list,color='mediumseagreen',linewidth=1.2,linestyle='dashed')
F.add_label(305.18,0.35,'N',color='seagreen',family='serif',weight='bold',size=14)
F.add_label(305.225,-0.13,'S',color='seagreen',family='serif',weight='bold',size=14)

F.show_markers(xhii,yhii,c='m',marker=hii_marker,s=70)
F.show_markers(xuchii,yuchii,c='m',marker=uchii_marker,s=70)
F.show_markers(xbrc,ybrc,c='m',marker=brc_marker,s=70)

F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/ncol_total_BIG.eps')


plt.rcParams.update({'font.size':14})
fig2 = plt.figure(figsize=(8,8))
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305-excitation.fits',figure=fig2)
F.set_theme('publication')
F.show_colorscale(cmap='jet',vmax=1.0)
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_font(size='medium', weight='medium', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')
F.colorbar.set_location('right')
F.colorbar.set_width(0.1)
F.colorbar.set_pad(0.07)
#F.colorbar.set_axis_label_text(r"$^{13}$CO Column Density (cm$^{-2}$)")
F.colorbar.set_axis_label_text(r"$^{13}$CO (3-2)/(2-1)")
F.colorbar.set_axis_label_pad(10)
F.colorbar.set_axis_label_rotation(0)
#F.set_title(r"$^{13}$CO (3-2)[LAsMA]/(2-1)[SEDIGISM]")
#F.axis_labels.hide()
#line_list = [np.array([[305.341,305.0688],[0.1210355,0.4414058]]),np.array([[305.3177,305.2177],[0.009924373,-0.1789645]])]
line_list = [np.array([[305.351,305.0788],[0.1210355,0.4414058]]),np.array([[305.3177,305.2177],[0.009924373,-0.1789645]])]
F.show_lines(line_list,color='k',linewidth=1.2,linestyle='dashed')
F.add_label(305.18,0.36,'N',color='k',family='serif',weight='bold',size=14)
F.add_label(305.225,-0.13,'S',color='k',family='serif',weight='bold',size=14)

#
# Mark Stars
#
import matplotlib.markers as markers

xstar = np.array([305.132,305.166,305.21,305.22,305.236,305.259,305.27,305.302,\
                 305.33,305.332,305.336,305.342,305.354,305.345,305.36,305.391,\
                 305.406,305.421,305.422,305.443,305.411,305.521,305.751,305.848,\
                 305.925])
ystar = np.array([0.0641604,-0.162111,0.0565067,0.104962,0.022534,0.226514,\
                 -0.0067114,0.0522876,0.104147,0.0289869,0.0539986,0.0777438,\
                 0.196546,0.0399662,0.11355,0.0383904,0.060785,0.10918,0.0148179,\
                 0.0126672,-0.822342,-0.826495,0.0308182,-0.0678552,-0.142517])
star_marker = markers.MarkerStyle(marker=(5,1,0),fillstyle="none")
F.show_markers(xstar,ystar,c='firebrick',marker=star_marker,s=30)

x_Danks = np.array([305.3384,305.3934])
y_Danks = np.array([+00.0719,+00.0874])
radii_Danks = np.array([0.018,0.026])
F.show_circles(x_Danks,y_Danks,radius=radii_Danks,edgecolor='k',facecolor='none',linewidths=1)

x_wr48 = 305.361
y_wr48 = +00.056
F.show_markers(x_wr48,y_wr48,c='k',marker=star_marker,s=60)

#
# Draw Beam and scalebar
#
from astropy import units as u
plt.rcParams.update({'font.size':12})
F.add_scalebar(0.1,label='6.6 pc',frame=True,linewidth=1)
F.add_beam()
F.beam.set_major(19*u.arcsecond)  # degrees
F.beam.set_minor(19*u.arcsecond)  # degrees
F.beam.set_edgecolor('black')
F.beam.set_facecolor('red')
#F.beam.set_linewidth(2)  # points
F.beam.set_corner('bottom right')

F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305-excitation-small.eps')
