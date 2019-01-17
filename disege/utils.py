# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:48:55 2019

@author: iago.costa
"""

from __future__ import division, print_function
import matplotlib.pyplot as plt
from fatiando.gravmag import prism, transform, euler
from fatiando.mesher import Prism
from fatiando import gridder, utils
import numpy as np
from scipy.interpolate import griddata
import matplotlib.patches as patches
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LightSource


    
def grid_data(x,y,z,area = None,cel = 50,method = 'cubic'):
    # Grid data for plot
    
    if not area:
        area = (x.min(),x.max(),y.min(),y.max())
    
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    print(numcols,numrows)
    x_y = np.c_[x.ravel(),y.ravel()]
    xi = np.linspace(area[0], area[1], numcols)
    yi = np.linspace(area[3], area[2], numrows)
    xi, yi = np.meshgrid(xi, yi)
    
    zi = griddata(x_y,z.ravel(),(xi,yi),method = method)
    
    return xi,yi,zi

def syn_data(inc, dec,syn_bounds,area,depth=-300,magnetization = 5,cel=200):
    # Build the synthetic Anomaly
    mag = utils.ang2vec(magnetization, inc, dec)
    model = []
    for i in range(len(syn_bounds)):
        model.append(Prism(syn_bounds[i][0], syn_bounds[i][1], syn_bounds[i][2], syn_bounds[i][3], syn_bounds[i][4], syn_bounds[i][5], {'magnetization': mag}))
    
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    x, y, z = gridder.regular(area[:4], (numcols, numrows), z=depth)

    #Calculate the Total Magnetic Field in nT
    mag = prism.tf(x, y, z, model, inc, dec)
    
    return x,y,mag,model


def syn_plot(xi,yi,z,area,syn_bounds,model,label,cel=10):
    #
    #Plotting synthetic data from xi,yi,z obtained from grid_data
    print(label)
    print(area)
    print(syn_bounds)
    fig = plt.figure(figsize=(20, 8))
    ax1 = fig.add_subplot('121')

    
    if model:
        for i in range(len(model)):
            rect = patches.Rectangle((model[i].get_bounds()[0],model[i].get_bounds()[2]),model[i].get_bounds()[1]-model[i].get_bounds()[0],model[i].get_bounds()[3]-model[i].get_bounds()[2],linewidth=1,edgecolor='black',facecolor='none')
            ax1.add_patch(rect)
            
    pcm = plt.contourf(xi, yi, z,50,vmin=z.min(), vmax=z.max(),cmap=plt.cm.rainbow,)
    #if you want to normalize in log > e.g. norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=z.min(), vmax=z.max())
    cbar = plt.colorbar(pcm,orientation='horizontal',fraction=0.046, pad=0.1)
    cbar.set_label('nT',fontsize = 12)
    ax1.set_xlabel('x (m)',fontsize = 12)
    ax1.set_ylabel('y (m)',fontsize = 12)
    ax1.set_title(label,fontsize = 16,style='italic')
    ax1.axis('scaled')
    
       
    # Plotting the 3D anomaly
    ax2 = fig.add_subplot(122,projection='3d')
    ax2.set_title('Synthetic Anomaly',fontsize = 16, style='italic')
    
    for i in range(len(model)):
        cube_definition = (model[i].get_bounds()[0],model[i].get_bounds()[2],model[i].get_bounds()[4]), (model[i].get_bounds()[0],model[i].get_bounds()[3],model[i].get_bounds()[4]), (model[i].get_bounds()[1],model[i].get_bounds()[2],model[i].get_bounds()[4]), (model[i].get_bounds()[0],model[i].get_bounds()[2],model[i].get_bounds()[5])
        edges,points = plot_cube(cube_definition)
    
    
        #cube_definition = (min_x,min_y,min_z), (min_x,max_y,min_z), (max_x,min_y,min_z), (min_x,min_y,max_z)
        faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
        faces.set_facecolor((0,0,1,0.1))
    
        ax2.add_collection3d(faces)
    
    # Plot the points themselves to force the scaling of the axes
    ax2.scatter(points[:,0], points[:,1], points[:,2], s=0)
    ax2.set_xlim(area[0],area[1])
    ax2.set_ylim(area[2],area[3])
    ax2.set_zlim(area[4],area[5])
    ax2.invert_zaxis()
    ax2.set_xlabel('x(m)')
    ax2.set_ylabel('y(m)')
    ax2.set_zlabel('z(m)')
    ax2.set_aspect('auto')

    plt.xlim(area[0],area[1])
    plt.ylim(area[2],area[3])
    plt.tight_layout()
    
    
def plot_cube(cube_definition):
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]

    return edges,points



def save_data(x,y,z,fname):
    f= open(fname,"w")
    for i in range(len(x)):
        f.write('%i\t%i\t%f\t\n' % (x[i],y[i],z[i]))
    f.close()


def load_xyz(fname, usecols = (0,1,3)):
    #Load xyz file from geosoft
    i = 0
    skip_list = []
    f= open(fname,"r")
    for line in f:
        if line.split()[0].startswith('Line'):
            skip_list.append(i)
        i = i+1
    
    x,y,z = np.loadtxt(fname,comments = ' Line',usecols = usecols,unpack = True)
    
    return x,y,z
    
    

def euler_deconv(inc, dec,syn_bounds,area,depth=-300,magnetization = 0.5,si=1.0,size = (1000, 1000),windows=(10, 10),proc_data='None'):
    
    # IN PROGRESSING #####
    
    
    mag = utils.ang2vec(magnetization, inc, dec)
    model = []
    for i in range(len(syn_bounds)):
        model.append(Prism(syn_bounds[i][0], syn_bounds[i][1], syn_bounds[i][2], syn_bounds[i][3], syn_bounds[i][4], syn_bounds[i][5], {'magnetization': mag}))
    #model = [Prism(syn_bounds[0]-1000, syn_bounds[1]-1000, syn_bounds[2]-1000, syn_bounds[3]-1000, syn_bounds[4], syn_bounds[5], {'magnetization': mag}),Prism(syn_bounds[0]+1000, syn_bounds[1]+1000, syn_bounds[2]+1000, syn_bounds[3]+1000, syn_bounds[4], syn_bounds[5], {'magnetization': mag})]
    
    
    cel = 200
    numcols, numrows = int((area[1] - area[0])/cel), int((area[3] - area[2])/cel)
    x, y, z = gridder.regular(area[:4], (numcols, numrows), z=depth)
    
    
    mag = prism.tf(x, y, z, model, inc, dec)
    #derivatives
    deriv_x = transform.derivx(x, y, mag, (numcols, numrows))
    deriv_y = transform.derivy(x, y, mag, (numcols, numrows))
    deriv_z = transform.derivz(x, y, mag, (numcols, numrows))
    #label = 'Total Magnetic Intensity (nT)'
    
    solver = euler.EulerDeconvMW(x, y, z, mag, deriv_x, deriv_y, deriv_z,
                             structural_index=si, windows=windows,
                             size=size)
    
    
    solver_expand = euler.EulerDeconvEW(x, y, z, mag, deriv_x, deriv_y, deriv_z,
                             structural_index=si, center = (-1250,0),sizes=np.linspace(300, 7000, 20))
    
    solver.fit()
    solver_expand.fit()
    #print('Kept Euler solutions after the moving window scheme:')
    #print(solver_expand.estimate_)
    return solver,solver_expand








def plot_data(xi,yi,z,label=None,shaded=False):
    #
    #Plotting data from xyz file
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot('111')
    
    if shaded:
        ls = LightSource(azdeg=315, altdeg=45)
        ve = 0.1
        #rgb = ls.shade(z,cmap=plt.cm.rainbow, blend_mode = 'overlay', vert_exag=1000, fraction = 0.3, dx = 200, dy = 200)
        plt.imshow(ls.hillshade(z,vert_exag=ve),cmap='gray')
        rgb = ls.shade(z,cmap=plt.cm.rainbow, blend_mode = 'hsv',vert_exag=ve)
        #pcm = plt.imshow(z,cmap=plt.cm.rainbow)
        plt.imshow(rgb,extent=[xi.min(),xi.max(),yi.min(),yi.max()])
        
        #if you want to normalize in log > e.g. norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=z.min(), vmax=z.max())
        
    else:
        #pcm = plt.contourf(xi, yi, z,50,vmin=z.min(), vmax=z.max(),cmap=plt.cm.rainbow)
        pcm = plt.imshow(z,extent=[xi.min(),xi.max(),yi.min(),yi.max()],cmap=plt.cm.rainbow,vmin = -300, vmax = 300)
    cbar = plt.colorbar(pcm,orientation='horizontal',fraction=0.046, pad=0.1)    
    cbar.set_label('nT',fontsize = 12)
    ax1.set_xlabel('x (m)',fontsize = 12)
    ax1.set_ylabel('y (m)',fontsize = 12)
    ax1.set_title(label,fontsize = 16,style='italic')
    plt.axis('scaled')
    
    
    plt.xlim(xi.min(),xi.max())
    plt.ylim(yi.min(),yi.max())
    plt.tight_layout()