#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:30:39 2021

@author: pmazumdar

Task: uses a modified version of the model from Watkins+ 2019\
      to calculate the radius of ionizaiton for a spherically
      expanding shell.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

save = True # set this to False if you dont want to save the output

def L(N,ne,nh,T):    # Eqn.23 from Watkins 2019 modified for spherical case
    a = 2.7e-13
    f = 4*np.pi*a*ne*ne
    f1 = 3*N/f
    t = T*31557600000000.0
    return (f1*(1-np.exp(-ne*ne*a*t/nh)))**(1/3)*u.cm.to(u.parsec) # returns radius in parsec
  
N_Danks2 = 1.03e50 # Calculated from list of stars in Danks 2 from Davies and Luminosity from Panagia paper
N_Danks1 = 5.586596185289926e+50
N_O6 = 3.1622776601683793e+49
N_WN = 1.2589254117941714e+49
N_WC = 1.0e49
N_other = 4*N_WC + 3*N_WN
nh=1e4 # upper limit based on Hindson 10 and 13
nh_alt = 1.5e3
ne = [60,100,500,5000] # electron densities

x = np.logspace(-5,1,100)   # time axis

colors = ['crimson','goldenrod','black','darkcyan','blue']

plt.rcParams.update({'font.size':11})
fig,ax = plt.subplots(figsize=(5.5,5))

for i,nex in enumerate(ne):
    #ax.plot(x,2*L(N_Danks2,nex,nh,x),color=colors[i],linewidth=0.7,linestyle='dashed')
    ax.plot(x,2*L(N_Danks1+N_Danks2+N_other,nex,nh,x),color=colors[i],linewidth=0.7,label='n$_e$ = '+str(nex)+'cm$^{-3}$')    
    #ax.plot(x,2*L(N_O6,nex,nh,x),color=colors[i],linewidth=0.7,linestyle='dashed',label='n$_e$ = '+str(nex)+'cm$^{-3}$')    
    ax.set_xscale('log')
    ax.set_yscale('log')
#ax.axvline(x=2,color='k',linewidth=0.7,linestyle='dotted')
ax.axvline(x=1,color='k',linewidth=0.7,linestyle='dotted')
#ax.text(2.0,0.7,r"Min. Age$_{\rm{Danks2}}$",rotation=90)
ax.text(1.0,6,r"Min. Age$_{   \rm{\.Danks1}}$",rotation=90)
#ax.text(1e-5,4.1,r"n$_{\rm{H}}$ = 10$^4$ cm$^{-3}$")
ax.set_ylabel("Ionized Diameter [pc]",fontsize=12)
ax.set_xlabel("Time [Myr]",fontsize=12)
#plt.legend(fancybox=True)
plt.legend(frameon=False)
plt.tight_layout()
if save:
  plt.savefig('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/ionized_length_time_evolution.eps',overwrite=True)