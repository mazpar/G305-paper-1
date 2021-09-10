# -*- coding: utf-8 -*-
"""
Script for plotting Integrated Maps

"""

import aplpy
import matplotlib.pyplot as plt
import numpy as np

#
# Moment-0 Maps
#
plt.rcParams.update({'font.size':16})
fig = plt.figure(figsize=(12,12))

# Load file
#F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_0.fits',figure=fig)
F = aplpy.FITSFigure('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO_mom_0.fits',figure=fig)

# Plot moment map
F.set_theme('publication')
F.show_colorscale(cmap='plasma',vmin=0)
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_font(size='medium', weight='medium', \
                    stretch='normal', family='sans-serif', \
                    style='normal', variant='normal')
F.colorbar.set_location('top')
F.colorbar.set_width(0.14)
F.colorbar.set_pad(0.07)
F.colorbar.set_axis_label_text(r"$^{13}$CO Moment-0 Map (K.km.s$^{-1}$)")
F.colorbar.set_axis_label_pad(30)
F.colorbar.set_axis_label_rotation(0)


# Mark Stars
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
F.show_markers(xstar,ystar,c='limegreen',marker=star_marker,s=80)

x_Danks = np.array([305.3384,305.3934])
y_Danks = np.array([+00.0719,+00.0874])
radii_Danks = np.array([0.018,0.026])
F.show_circles(x_Danks,y_Danks,radius=radii_Danks,edgecolor='red',facecolor='none',linewidths=2)

x_wr48 = 305.361
y_wr48 = +00.056
F.show_markers(x_wr48,y_wr48,c='red',marker=star_marker,s=120)


# Draw Beam and scalebar
from astropy import units as u
F.add_scalebar(0.1,label='6.6 pc',frame=True,linewidth=2)
F.add_beam()
F.beam.set_major(19*u.arcsecond)  # degrees
F.beam.set_minor(19*u.arcsecond)  # degrees
F.beam.set_edgecolor('black')
F.beam.set_facecolor('red')
#F.beam.set_linewidth(2)  # points
F.beam.set_corner('bottom right')


# Draw rectangles marking distinct regions
x_c = [305.87,305.62,305.5,305.55,305.38,305.21,305.10,305.18,305.53] # centre x-coordinate
y_c = [-0.13,0.13,-0.06,0.38,0.24,0.3,0.13,-0.08,-0.35] # centre y-coordinate
width = [0.25,0.17,0.29,0.2,0.15,0.18,0.12,0.25,0.45]
height = [0.35,0.2,0.17,0.2,0.22,0.25,0.15,0.15,0.3]
angle = [0.0,45.0,30.0,0.0,0.0,0.0,0.0,25.,0.]
shift13 = [0.1,-0.1,-0.1,-0.08,0,0,0,0,0.1]
layer = ['I','II','III','IV','V','VI','VII','VIII','IX'] # layer names

# pick colors for the rectangles.
#colors = ['white','white','white','white','white','white','white','white','white'] # 12CO
colors = ['k','k','k','k','k','k','k','k','k'] # 13CO

F.show_rectangles(x_c, y_c, width, height, angle=angle, zorder=None, coords_frame='world',color=colors,linewidth=1.5,linestyle='dashed')


#for i, layer_name in enumerate(layer):   # 12CO
#	F.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 

for i, layer_name in enumerate(layer):   # 13CO
	F.add_label(x_c[i]-shift13[i],y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 

#F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_12CO_mom_0_BIG.eps') # 12CO
F.save('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_13CO_mom_0_BIG.eps') # 13CO
