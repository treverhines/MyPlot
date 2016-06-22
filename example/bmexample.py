#!/usr/bin/env python
import numpy as np
from  myplot.basemap import Basemap
import matplotlib.pyplot as plt

lon_0 = -116.0
lat_0 = 33.4

llcrnrlon = -118.0
urcrnrlon = -114.0
llcrnrlat = 31.0
urcrnrlat = 35.0

fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
bm = Basemap(lon_0=lon_0,lat_0=lat_0,
             llcrnrlon=llcrnrlon,
             llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon,
             urcrnrlat=urcrnrlat,
             projection='tmerc',resolution='h',ax=ax)

#ext = ax.get_window_extent()
#def get_axis_rect(ax,fig):
#  ext = ax.get_window_extent()
#  left = (ext.xmin/fig.dpi)/fig.get_figwidth()
#  bottom = (ext.ymin/fig.dpi)/fig.get_figheight()
#  width = ((ext.xmax-ext.xmin)/fig.dpi)/fig.get_figwidth()
#  height = ((ext.ymax-ext.ymin)/fig.dpi)/fig.get_figheight()
#  return left,bottom,width,height

#def make_inset_ax(ax,fig,):
#  rect = get_axis_rect(ax,fig)  

  
#print(get_axis_rect(ax,fig))
#bm.drawcoastlines()
#bm.etopo(scale=4.0)
lonmin = bm.lonmin
lonmax = bm.lonmax
# add a 10% buffer
lonmin -= 0.1*(lonmax - lonmin)
lonmax += 0.1*(lonmax - lonmin)

latmin = bm.latmin
latmax = bm.latmax
latmin -= 0.1*(latmax - latmin)
latmax += 0.1*(latmax - latmin)

lons = np.linspace(lonmin,lonmax,100)
lats = np.linspace(latmin,latmax,100)
lons,lats = np.meshgrid(lons,lats)
lonlats = np.array([lons.flatten(),lats.flatten()]).T
data = np.sin(lons).flatten()

bm.drawscalar(data,lonlats,resolution=1000)
#bm.drawtopography(resolution=1000.0)
bm.drawminimap([0.1,0.1,0.2,0.2])
plt.show()

