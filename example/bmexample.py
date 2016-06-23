#!/usr/bin/env python
import numpy as np
from  myplot.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm

lon_0 = -116.0
lat_0 = 33.4
llcrnrlon = -118.0
urcrnrlon = -114.0
llcrnrlat = 31.0
urcrnrlat = 35.0

bm = Basemap(lon_0=lon_0,lat_0=lat_0,
             llcrnrlon=llcrnrlon,
             llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon,
             urcrnrlat=urcrnrlat,
             projection='tmerc',resolution='i')


fig,ax = plt.subplots()
x = np.random.random((1000,2))
x = bm.axes_to_geodetic(x,ax)
val = np.sin(2*x[:,0])*np.cos(2*x[:,1])

bm.drawscalar(val,x,resolution=500,zorder=0,topography=True,ax=ax,cmap=matplotlib.cm.viridis)
bm.drawtopography(resolution=500,alpha=0.2,ax=ax)
bm.drawcoastlines()

fig,ax = plt.subplots()
bm.drawscalar(val,x,resolution=500,zorder=0,topography=True,ax=ax,cmap=matplotlib.cm.viridis)
bm.drawtopography(resolution=500,alpha=0.0,ax=ax)
bm.drawcoastlines()

fig,ax = plt.subplots()
bm.drawtopography(resolution=500,alpha=0.2,ax=ax,cmap=matplotlib.cm.viridis)
bm.drawcoastlines()
plt.show()

