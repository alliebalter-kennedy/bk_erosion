# greenland_plot_erosion_result.py - BK Erosion test - Trying to make this commit into a new branch. and now allie is trying to add something at the same time as benjamin
import netCDF4 as nc
#import matplotlib.pyplot as plt
import numpy as np
import matlab
import matlab.engine as mlab
from matplotlib import pyplot as plt, cm
from matplotlib import colors
# Benjamin added this line.

fn='/Users/keislingldeo/Documents/lamont21-22/greenland_erosion/data/rip33_gdgt_1.75_0.5.nc'


#al=np.loadtxt("./rip33_gdgt_1.75_0.5/aldepthsum_final.txt")
#be=np.loadtxt("./rip33_gdgt_1.75_0.5/bedepthsum_final.txt")

data=np.loadtxt("./rip33_gdgt_1.75_0.5/cdepthsum_final.txt")


##PLOT THE TOTAL EROSION 
#plt.imshow(erosionsum, origin='lower', cmap=cm.viridis, norm=colors.LogNorm())
#plt.clim([1e-1,1e3])
#plt.colorbar()
#plt.show()



##PLOT THE BERYLLIUM RESULTS
plt.imshow(data, origin='lower', cmap=cm.rainbow)
#plt.imshow(data, origin='lower', cmap=cm.rainbow, norm=colors.LogNorm())
#plt.imshow(al/be, origin='lower', cmap=cm.rainbow)
plt.clim([0,1e3])
plt.colorbar()
plt.show()


#fig, (ax1,ax2) = plt.subplots(2,3)
#ax1.imshow(data1, origin='lower', cmap=cm.rainbow)
#ax2.imshow(data2, origin='lower', cmap=cm.rainbow)

#plt.show()
