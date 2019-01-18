# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:48:55 2019

@author: iago.costa
"""

from __future__ import division, print_function
import matplotlib.pyplot as plt
#DISEGE library
from disege import processing,utils



##### Create some synthetic magnetic data #########
    
def main():
    # Inclination and Declination ##
    inc, dec = 90, 0 # B = 57.000 - parameters at the pole 
    #inc, dec = -9.5,-20.4  # B = 24.768 - parameters at Amazonian terrain
    
    # Synthetic anomalie boundaries: min_x,max_x,min_y,max_y,min_z,max_z,magnetization ##
    #syn_bounds = [(-1500,-1000,-2500,2500,100,1100),(0,500,-2500,2500,300,1300),(1500,2000,-2500,2500,500,1500)]   # Example of three retangles
    syn_bounds = [(-1000,1000,-1000,1000,100,1100)]
    
    # Model Area ##
    area = (-2000,2000,-2000,2000,0,2000)
    
    # Build the Synthetic Anomaly #
    x,y,mag,model = utils.syn_data(inc, dec,syn_bounds,area,magnetization = 5)
    
    # Save Data into a xyz file to import in Geosoft - Opitionally
    utils.save_data(x,y,mag,fname = 'mag_disege_pole.xyz')
    
    #Processing data - Optionally
    asa = processing.gdo_mag(x,y,mag,area,theta=0,phi=45)
    
    # Grid data
    xi,yi,z = utils.grid_data(x,y,asa,area)
    
    #Plotting synthetic data
    utils.syn_plot(xi,yi,z,area,syn_bounds,model,label='Analytic Signal (nT/m)',clabel = 'teste')
    
    #Show the final plot
    plt.show()
    
main()   

    
    





















