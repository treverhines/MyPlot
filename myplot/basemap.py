#!/usr/bin/env python
from mpl_toolkits.basemap import Basemap as _Basemap
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm
from pegtop import set_shade
from pegtop import hillshade
from rbf.interpolant import RBFInterpolant
import rbf.basis
import matplotlib.patches as patches

class Basemap(_Basemap):
  def drawtopography(self,resolution=200,cmap=matplotlib.cm.gray,vmin=None,vmax=None,**kwargs):
    try:
      etopodata = h5py.File('/cmld/data5/hinest/Projects/MyPlot/data/ETOPO1.h5')
    except IOError:
      print('cannot open ETOPO1')
      return

    lons = etopodata['longitude'][...]
    lats = etopodata['latitude'][...]

    lonmin = self.lonmin   
    lonmax = self.lonmax
    # add a 10% buffer
    lonmin -= 0.1*(lonmax - lonmin)
    lonmax += 0.1*(lonmax - lonmin)

    latmin = self.latmin
    latmax = self.latmax
    latmin -= 0.1*(latmax - latmin)
    latmax += 0.1*(latmax - latmin)


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
    etopodata.close() 
    nx = resolution#int((self.xmax-self.xmin)/resolution)+1
    ny = resolution#int((self.ymax-self.ymin)/resolution)+1
    topodat = self.transform_scalar(data,lons,lats,nx,ny,order=3,masked=True)
    rgb = set_shade(topodat,cmap=cmap,scale=10.0,azdeg=300.0,altdeg=45.0,vmin=vmin,vmax=vmax)
    self.imshow(rgb,**kwargs)

  def drawscalar(self,val,lonlat,topography=False,resolution=200,penalty=0.0,cmap=matplotlib.cm.seismic,**kwargs):
    lonmin = self.lonmin   
    lonmax = self.lonmax
    # add a 10% buffer
    lonmin -= 0.1*(lonmax - lonmin)
    lonmax += 0.1*(lonmax - lonmin)

    latmin = self.latmin
    latmax = self.latmax
    latmin -= 0.1*(latmax - latmin)
    latmax += 0.1*(latmax - latmin)

    # form data interpolant        
    A = RBFInterpolant(lonlat,val,order=0,basis=rbf.basis.phs1,penalty=penalty)
    nx = resolution#int((self.xmax-self.xmin)/resolution)
    ny = resolution#int((self.ymax-self.ymin)/resolution)
  
    lons = np.linspace(lonmin,lonmax,nx)
    lats = np.linspace(latmin,latmax,ny)
    latitp,lonitp = np.meshgrid(lats,lons)
    latitp = latitp.flatten()
    lonitp = lonitp.flatten()
    positp = np.array([lonitp,latitp]).T
    valitp = A(positp)
    valitp = np.reshape(valitp,(nx,ny)).T
    data = self.transform_scalar(valitp,lons,lats,nx,ny,order=3)

    # form topography based intensity if topography is True
    if topography:
      try:
        etopodata = h5py.File('/cmld/data5/hinest/Projects/MyPlot/data/ETOPO1.h5')
      except IOError:
        print('cannot open ETOPO1')
        return

      lons = etopodata['longitude'][...]
      lats = etopodata['latitude'][...]
    
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
      topo_data = etopodata['elevation'][ystart:yend,xstart:xend]
      lonsgrid,latsgrid = np.meshgrid(lons,lats)
      topo_data = self.transform_scalar(topo_data,lons,lats,nx,ny,order=3,masked=True)
      intensity = hillshade(topo_data,scale=10.0,azdeg=300.0,altdeg=45.0)
      etopodata.close() 
      rgb = set_shade(data,intensity=intensity,cmap=cmap)
      return self.imshow(rgb,cmap=cmap,**kwargs)
  


    else:    
      return self.imshow(data,cmap=cmap,**kwargs)
  
  def grid(self,spacing=1.0,**kwargs):
    self.drawmeridians(np.arange(np.floor(self.llcrnrlon),
                     np.ceil(self.urcrnrlon),spacing),
                     labels=[0,0,0,1],
                     **kwargs)
    self.drawparallels(np.arange(np.floor(self.llcrnrlat),
                       np.ceil(self.urcrnrlat),spacing),
                       labels=[1,0,0,0],
                       **kwargs)
    return

  def axes_to_geodetic(self,pnt,ax):
    pnt_disp = ax.transAxes.transform(pnt)
    pnt_data = ax.transData.inverted().transform(pnt_disp)
    pnt_geo = self(pnt_data[0],pnt_data[1],inverse=True)
    return np.array(pnt_geo)

  def geodetic_to_axes(self,pnt,ax):
    pnt_data = self(pnt[0],pnt[1])
    pnt_disp = ax.transData.transform(pnt_data)
    pnt_ax = ax.transAxes.inverted().transform(pnt_disp)
    return pnt_ax
      
  def drawminimap(self,rect,ax=None,fig=None,scale=5.0):
    if ax is None:
      ax = plt.gca()
    if fig is None:
      fig = ax.get_figure()

    miniax = fig.add_axes(rect)
    width = scale*(self.xmax - self.xmin)
    height = scale*(self.ymax - self.ymin)
    lon_0 = self.lonmin + (self.lonmax - self.lonmin)/2.0
    lat_0 = self.latmin + (self.latmax - self.latmin)/2.0
    bm = Basemap(lon_0=lon_0,lat_0=lat_0,width=width,height=height,
                 ax=miniax,projection=self.projection,resolution='i')
    bm.drawcoastlines()
    bm.drawcountries()
    bm.drawstates()
    rect_min,rect_max = bm(self.llcrnrlon,self.llcrnrlat)
    rect_width = (self.xmax - self.xmin)
    rect_height = (self.ymax - self.ymin)
    miniax.add_patch(patches.Rectangle((rect_min,rect_max),
                                        rect_width,rect_height,
                                        facecolor='red',alpha=0.4,
                                        edgecolor='none'))
    plt.sca(ax)
    return
  
