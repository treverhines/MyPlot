#!/usr/bin/env python
import scipy.interpolate
import numpy as np
import myplot.cm
from myplot.xsection import cmap_to_rgba
import h5py
import mayavi.mlab 

this_dir = '/'.join(__file__.split('/')[:-1])

def draw_topography(basemap,zscale=1,cmap=myplot.cm.etopo1,resolution=0.02,**kwargs):
  #etopodata = h5py.File('/cmld/data5/hinest/PyMods/myplot/data/ETOPO1.h5')
  try:
    #etopodata = h5py.File('%s/data/ETOPO1.h5' % this_dir)
    etopodata = h5py.File('/cmld/data5/hinest/Projects/MyPlot/data/ETOPO1.h5')
  except IOError:
    print('cannot load ETOPO1.h5')
    return
  
  lons = etopodata['longitude'][...]
  lats = etopodata['latitude'][...]

  lonmin = basemap.lonmin
  lonmax = basemap.lonmax
  latmin = basemap.latmin
  latmax = basemap.latmax

  xidx = (lons > lonmin) & (lons < lonmax)
  xidx = xidx.nonzero()[0]
  xstart = xidx[0]
  xend = xidx[-1]

  yidx = (lats > latmin) & (lats < latmax)
  yidx = yidx.nonzero()[0]
  ystart = yidx[0]
  yend = yidx[-1]

  lons = lons[xstart:xend]
  lats = lats[ystart:yend]
  data = etopodata['elevation'][ystart:yend,xstart:xend]

  #data = topoin.data                                                                                              
  itp = scipy.interpolate.interp2d(lons,lats,data,kind='cubic')

  lonitp = np.arange(min(lons),max(lons),resolution)
  latitp = np.arange(min(lats),max(lats),resolution)
   
  dataitp = itp(lonitp,latitp)
  longrid,latgrid = np.meshgrid(lonitp,latitp)
  xgrid,ygrid = basemap(longrid,latgrid)
  vmin = -8210.0*zscale
  vmax = 7000.0*zscale
  m = mayavi.mlab.mesh(xgrid,ygrid,zscale*dataitp,vmin=vmin,vmax=vmax,
                       **kwargs)
  rgba = cmap_to_rgba(cmap)
  m.module_manager.scalar_lut_manager.lut.table = rgba
  mayavi.mlab.draw()



