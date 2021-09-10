#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:03:59 2020

@author: pmazumdar

Task: Plotting Channel Maps for a data_cube
"""
################# Plotting Average Spectra over a masked cube #################
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
import matplotlib
import aplpy
from spectral_cube import SpectralCube


cube = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_resample.fits')
################ Initialize Figure ########################
fig = plt.figure(figsize=(30, 15))
matplotlib.rcParams.update({'font.size': 18})



F1 = aplpy.FITSFigure(cube[60,:,:].hdu, figure=fig,subplot=[0.1,0.1,0.2,0.4] )
F1.set_theme('publication')
F1.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F1.add_colorbar()
F1.colorbar.set_axis_label_text("[K]")
F2 = aplpy.FITSFigure(cube[70,:,:].hdu, figure=fig,subplot=[0.35,0.1,0.2,0.4] )
F2.set_theme('publication')
F2.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F2.add_colorbar()
F2.colorbar.set_axis_label_text("[K]")
F3 = aplpy.FITSFigure(cube[80,:,:].hdu, figure=fig,subplot=[0.6,0.1,0.2,0.4] )
F3.set_theme('publication')
F3.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F3.add_colorbar()
F3.colorbar.set_axis_label_text("[K]")
F4 = aplpy.FITSFigure(cube[110,:,:].hdu, figure=fig,subplot=[0.1,0.55,0.2,0.4] )
F4.set_theme('publication')
F4.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F4.add_colorbar()
F4.colorbar.set_axis_label_text("[K]")
F5 = aplpy.FITSFigure(cube[100,:,:].hdu, figure=fig,subplot=[0.35,0.55,0.2,0.4] )
F5.set_theme('publication')
F5.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F5.add_colorbar()
F5.colorbar.set_axis_label_text("[K]")
F6 = aplpy.FITSFigure(cube[90,:,:].hdu, figure=fig,subplot=[0.6,0.55,0.2,0.4] )
F6.set_theme('publication')
F6.show_colorscale(vmin=0,vmid=8,vmax=20,stretch='sqrt')
F6.add_colorbar()
F6.colorbar.set_axis_label_text("[K]")

F1.axis_labels.set_font(size='large')
F1.tick_labels.set_font(size='large')
F2.axis_labels.hide()
F2.tick_labels.hide()
F3.axis_labels.hide()
F3.tick_labels.hide()
F4.axis_labels.hide()
F4.tick_labels.hide()
F5.axis_labels.hide()
F5.tick_labels.hide()
F6.axis_labels.hide()
F6.tick_labels.hide()

matplotlib.rcParams.update({'font.size': 18})

F1.add_label(0.15, 0.9, '-25 km/s',color='white', relative=True)
F2.add_label(0.15, 0.9, '-30 km/s',color='white', relative=True)
F3.add_label(0.15, 0.9, '-35 km/s',color='white', relative=True)
F4.add_label(0.15, 0.9, '-50 km/s',color='white', relative=True)
F5.add_label(0.15, 0.9, '-45 km/s',color='white', relative=True)
F6.add_label(0.15, 0.9, '-40 km/s',color='white', relative=True)


#
# Draw rectangles
#

x_c = [305.87,305.62,305.5,305.55,305.39,305.21,305.10,305.18,305.52]
y_c = [-0.13,0.13,-0.06,0.38,0.24,0.3,0.13,-0.08,-0.35]
width = [0.25,0.17,0.29,0.2,0.13,0.18,0.12,0.25,0.45]
height = [0.35,0.2,0.17,0.2,0.2,0.25,0.15,0.15,0.3]
angle = [0.0,45.0,30.0,0.0,0.0,0.0,0.0,25.,0.]
layer = ['I','II','III','IV','V','VI','VII','VIII','IX']
colors = ['white','white','white','white','white','white','white','white','white']

F1.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
F2.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
F3.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
F4.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
F5.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
F6.show_rectangles(x_c, y_c, width, height, angle=angle,color=colors,linewidth=0.8,linestyle='dashed')
for i, layer_name in enumerate(layer):
	F1.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
	F2.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
	F3.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
	F4.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
	F5.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
	F6.add_label(x_c[i]+(angle[i]*width[i]*3.14/360),y_c[i],layer_name,color=colors[i],family='serif',weight='bold',size=14) 
fig.canvas.draw()
#fig.tight_layout()

fig.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/channel_maps_publication_new.eps',dvi=400)
