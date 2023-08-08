import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
import re
import math
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator

path = sys.argv[1]
band_name = sys.argv[2]



def get_data(fits_path):
    data  = fits.getdata(fits_path,ext=0)
    data  = data[::-1,:]
    nx,ny = data.shape
    
    hdulist = fits.open(fits_path)
    xmin  = float(hdulist[0].header['xmin'])
    xmax  = float(hdulist[0].header['xmax'])
    ymin  = float(hdulist[0].header['ymin'])
    ymax  = float(hdulist[0].header['ymax'])
    hdulist.close()

    return data,xmin,xmax,ymin,ymax,nx,ny



light_path = os.path.join(path,'output/'+band_name+'_lens_light_super.fits')
light,l_xmin,l_xmax,l_ymin,l_ymax,l_nx,l_ny = get_data(light_path)

mass_path = os.path.join(path,'output/lens_kappa_star_super.fits')
mass,m_xmin,m_xmax,m_ymin,m_ymax,m_nx,m_ny = get_data(mass_path)

# Interpolate mass on the light pixel coordinates
x_step = np.absolute(m_xmax - m_xmin)/m_nx
m_xx = [ (m_xmin + x_step/2.0 + x_step*i) for i in range(0,m_nx) ]
y_step = np.absolute(m_ymax - m_ymin)/m_ny
m_yy = [ (m_ymin + y_step/2.0 + y_step*i) for i in range(0,m_ny) ] 

x_step = np.absolute(l_xmax - l_xmin)/l_nx
l_xx = [ (l_xmin + x_step/2.0 + x_step*i) for i in range(0,l_nx) ]
y_step = np.absolute(l_ymax - l_ymin)/l_ny
l_yy = [ (l_ymin + y_step/2.0 + y_step*i) for i in range(0,l_ny) ] 

new_mass = mass

#interp = RegularGridInterpolator((m_xx,m_yy), mass)
#X,Y = np.meshgrid(l_xx,l_yy)
#new_mass = interp((X,Y))

# if l_xmin != m_xmin or l_xmax != m_xmax or l_ymin != m_ymin or l_ymax != m_ymax:
#     print("Ranges must be the same!")
#     sys.exit()
# if l_nx != m_nx or l_ny != m_ny:
#     print("Dimensions must be the same!")
#     sys.exit()




print(new_mass.shape)
print(light.shape)
diff = np.divide(light,new_mass)


    

mycmap_div    = 'Spectral'
mycmap_seq    = 'YlGnBu'
fig,axs = plt.subplots(1,3,figsize=(20,7))


im = axs[0].imshow(light,interpolation='none',extent=[l_xmin,l_xmax,l_ymin,l_ymax])
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right",size="5%",pad=0.05)
plt.colorbar(im,cax=cax,format="%6.4f")

#im = axs[1].imshow(mass,interpolation='none',extent=[m_xmin,m_xmax,m_ymin,m_ymax])
im = axs[1].imshow(mass,interpolation='none',extent=[l_xmin,l_xmax,l_ymin,l_ymax])
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right",size="5%",pad=0.05)
plt.colorbar(im,cax=cax,format="%6.4f")

im = axs[2].imshow(diff,interpolation='none',extent=[l_xmin,l_xmax,l_ymin,l_ymax],cmap=mycmap_div)
divider = make_axes_locatable(axs[2])
cax = divider.append_axes("right",size="5%",pad=0.05)
plt.colorbar(im,cax=cax,format="%6.4f")


# Zoom in
fac = 0.2
xmin = -fac
xmax = fac
ymin = -fac
ymax = fac
for ax in axs:
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)


plt.tight_layout()
plt.savefig("mass_to_light.png",bbox_inches='tight')
