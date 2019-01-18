# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:48:55 2019

@author: iago.costa
"""

from __future__ import division, print_function
from fatiando.gravmag import transform
import numpy as np


    
def thd_mag(x,y,mag,area,cel=200):
    # Total Horizontal Gradient (THD)
    #(Cordell and Grauch 1985)

    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    
    #applying the equations
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    
    return thd
    

def tdr_mag(x,y,mag,area,cel=200):
    # Tilt Angle Filter (TDR or ISA)
    #(Miller and Singh 1994)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    #applying the equations
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    tdr = np.arctan(deriv_z/thd)
    
    return tdr
    
def thd_tdr_mag(x,y,mag,area,cel=200):
    # Total Horizontal Derivative of the Tilt Angle (THD_TDR)
    #(Verduzco et al. 2004)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    #applying the equations
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    tdr = np.arctan(deriv_z/thd)
    derivx_tdr = transform.derivx(x, y, tdr, (numcols, numrows))
    derivy_tdr = transform.derivy(x, y, tdr, (numcols, numrows))
    
    thd_tdr = np.sqrt(derivx_tdr**2 + derivy_tdr**2)
    
    return thd_tdr
    
def asa_mag(x,y,mag,area,cel=200):
    # Analytical Signal or Total Gradient (ASA or GT)
    #(Nabighian 1972; Roest et al. 1992)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    #applying the equations
    asa = np.sqrt(deriv_x**2+deriv_y**2+deriv_z**2)
    
    return asa
    
def theta_mag(x,y,mag,area,cel=200):
    # Theta Map
    #(Wijns et al. 2005)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    #applying the equations
    asa = np.sqrt(deriv_x**2+deriv_y**2+deriv_z**2)
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    theta = (thd/asa)
    return theta

def tdx_mag(x,y,mag,area,cel=200):
    #TDX Filter
    #(Cooper and Cowan 2006)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    #applying the equations
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    tdx = np.arctan(thd/np.abs(deriv_z))
    
    return tdx
    
def tahg_mag(x,y,mag,area,cel=200):
    #TAHG Filter
    #(Ferreira et al. 2013)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
        
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    derivx_thd = transform.derivx(x,y,thd, (numcols,numrows))
    derivy_thd = transform.derivy(x,y,thd, (numcols,numrows))
    derivz_thd = transform.derivz(x,y,thd, (numcols,numrows))
    
    #applying the equations
    tahg = np.arctan(derivz_thd/np.sqrt(derivx_thd**2 + derivy_thd**2))
    
    return tahg
    
def stdr_mag(x,y,mag,area,cel=200,M=50000):
    #STDR Map
    #(Nasuti et al. 2018)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    deriv2_z = transform.derivz(x,y,deriv_z, (numcols,numrows))
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    derivx_thd = transform.derivx(x,y,thd, (numcols,numrows))
    derivy_thd = transform.derivy(x,y,thd, (numcols,numrows))
    
    #applying the equations
    stdr = np.arctan(M*deriv2_z/np.sqrt(derivx_thd**2+derivy_thd**2))
    
    return stdr
    

def thd_stdr_mag(x,y,mag,area,cel=200,M=50000):
    #THD_STDR Map
    #(Nasuti et al. 2018)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    
    deriv2_z = transform.derivz(x,y,deriv_z, (numcols,numrows))
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    derivx_thd = transform.derivx(x,y,thd, (numcols,numrows))
    derivy_thd = transform.derivy(x,y,thd, (numcols,numrows))
    
    #applying the equations
    stdr = np.arctan(M*deriv2_z/np.sqrt(derivx_thd**2+derivy_thd**2))
    derivx_stdr = transform.derivx(x,y,stdr,(numcols,numrows))
    derivy_stdr = transform.derivy(x,y,stdr,(numcols,numrows))
    thd_stdr = np.sqrt(derivx_stdr**2+derivy_stdr**2)
    
    return thd_stdr
    
def asta_mag(x,y,mag,area,cel=200):
    #Analytic Signal of Tilt Angle
    #(Ansari and Alamdar, 2011)
    
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))

    #Calculating TDR and its derivatives
    thd = np.sqrt(deriv_x**2 + deriv_y**2)
    tdr = np.arctan(deriv_z/thd)
    deriv_tdr_x = transform.derivx(x, y, tdr, (numcols, numrows))
    deriv_tdr_y = transform.derivy(x, y, tdr, (numcols, numrows))
    deriv_tdr_z = transform.derivz(x, y, tdr, (numcols, numrows))
    
    #applyng the equations
    asta = np.sqrt(deriv_tdr_x**2+deriv_tdr_y**2+deriv_tdr_z**2)
    
    return asta


def gdo_mag(x,y,mag,area,theta,phi,cel=200):
    #Generalised derivative operator
    #(Cooper and Cowan 2011)
    #Calculating the derivatives through Fatiando library
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    
    gdo = ((deriv_x*np.sin(theta) + deriv_y*np.cos(theta))*np.cos(phi) + deriv_z*np.sin(phi))/np.sqrt(deriv_x**2+deriv_y**2+deriv_z**2)
    
    return gdo
    
    