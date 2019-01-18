# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:48:55 2019

@author: iago.costa
"""

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
#DISEGE library
from disege import processing,utils



##### Plot Real Data from xyz file #########
    
def main():
    # Load xyz data
    x,y,z = utils.load_xyz(fname = 'MAG_cinz_reproj.xyz')
    
    #Grid Data
    xi,yi,zi = utils.grid_data(x,y,z,cel = 125)
    
    #Plotting data
    utils.plot_data(xi,yi,zi,label='teste',clabel='teste',shaded=False)
    
    #Show the final plot
    plt.show()
    
    
    #Plotting Real Data from .xyz file - UTM
    #xi,yi,z = grid_data(data_file='../MVI_xyz.xyz',area=area)
    #fig_plot(syn_bounds,area,xi,yi,z,label='Magnetic Vector Inversion (MVI)')
    
    

main()   

    
    





















