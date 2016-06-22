#!/usr/bin/env python
from misc import rotation3D
from mpl_toolkits.mplot3d import Axes3D as _Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib import cm
import numpy as np
import cm as tcm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from modest.tikhonov import Perturb

class Arrow3D(FancyArrowPatch):
  '''                                   
  creates an artist object for a 3D arrow      
  '''
  def __init__(self,xs,ys,zs,*args,**kwargs):
    '''                               
    arguments:             

      xs: the x coordinates of both sides of the arrow        

      ys: the y coordinates of both sides of the arrow      

      zs: the z coordinates of both sides of the arrow       

    '''
    FancyArrowPatch.__init__(self,(0,0),(0,0),*args,**kwargs)
    self.verts = xs,ys,zs

  def draw(self,renderer):
    xp,yp,zp = proj3d.proj_transform(*self.verts,M=renderer.M)
    self.set_positions((xp[0],yp[0]),(xp[1],yp[1]))
    FancyArrowPatch.draw(self,renderer)


class Axes3D(_Axes3D):
  def __init__(self,*args,**kwargs):
    if (len(args) > 1) | kwargs.has_key('rect'):
      _Axes3D.__init__(self,*args,**kwargs)

    else:
      rect = (0.1,0.1,0.8,0.8)
      _Axes3D.__init__(self,*args,rect=rect,**kwargs)

    self.cross_section_clim_set = False
    self.pbaspect = np.array([1.0,1.0,1.0])
    

  def cross_section(self,func,anchor,
                    zrot,yrot,xrot,
                    length,width,
                    Nl=20,Nw=20,
                    cmap=tcm.slip,
                    lw=1.0,
                    clim=None, 
                    func_args=None,
                    func_kwargs=None,
                    **kwargs):
    if func_args is None:
      func_args = ()

    if func_kwargs is None:
      func_kwargs = {}

    x = np.linspace(0,length,Nl)
    y = np.linspace(0,-width,Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    R = rotation3D(zrot,yrot,xrot)
    p = np.concatenate((x[None,:,:],
                        y[None,:,:],
                        z[None,:,:]),
                        axis=0)

    p = np.einsum('ij,jkl->kli',R,p)
    p += anchor
    p = np.reshape(p,(Nl*Nw,3))   
    c = func(p,*func_args,**func_kwargs)
    c = np.reshape(c,(Nl,Nw))
    p = np.reshape(p,(Nl,Nw,3))

    if not self.cross_section_clim_set:
      self.sm = cm.ScalarMappable(cmap=cmap)
      self.sm.set_array(c)
      if (clim is not None):
        self.sm.set_clim(clim[0],clim[1])

      self.cross_section_clim_set = True

    cnorm = self.sm.norm(c)#(c - clim[0])/(clim[1] - clim[0])
    s = self.plot_surface(p[:,:,0],
                          p[:,:,1],
                          p[:,:,2],
                          shade=False,
                          facecolors=cmap(cnorm),
                          rstride=1,cstride=1,
                          **kwargs)

    idx1 = np.array([[0,-1],[0,-1]])
    idx2 = np.array([[0,0],[-1,-1]])
    w = self.plot_wireframe(p[idx1,idx2,0],
                            p[idx1,idx2,1],
                            p[idx1,idx2,2],
                            color='k',lw=lw)

    return self.sm

  def vector_cross_section(self,func,anchor,
                           zrot,yrot,xrot,
                           length,width,
                           Nl=20,Nw=20,
                           lw=1.0,
                           minmag=0.05, 
                           arrow_length=1000.0,
                           func_args=None,
                           func_kwargs=None,
                           **kwargs):
    if func_args is None:
      func_args = ()

    if func_kwargs is None:
      func_kwargs = {}

    x = np.linspace(0,length,Nl)
    y = np.linspace(0,-width,Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    R = rotation3D(zrot,yrot,xrot)
    p = np.concatenate((x[None,:,:],
                        y[None,:,:],
                        z[None,:,:]),
                        axis=0)

    p = np.einsum('ij,jkl->kli',R,p)
    p += anchor
    p = np.reshape(p,(Nl*Nw,3))   
    c = func(p,*func_args,**func_kwargs)
    cmag = np.sqrt(c[:,0]**2 + c[:,1]**2 + c[:,2]**2)
    #c = np.ma.masked_array(c,mask=cmag<minmag)
    c *= arrow_length
    ends = p + c
    mscale = kwargs.pop('mutation_scale',1.0)
    for i in range(len(p)):
      if cmag[i] > minmag:
        artist = Arrow3D(*zip(p[i],ends[i]),mutation_scale=cmag[i]*mscale,**kwargs)
        self.add_artist(artist)


  def set_pbaspect(self,x):
    self.pbaspect = np.asarray(x)


  def draw_basis(self,
                 center=None,
                 length=1.0,
                 labels=None,
                 fontsize=16,
                 **kwargs):
    if center is None:
      center = np.array([0.0,0.0,0.0])

    center = np.asarray(center)
    if labels is None:
      labels = ['x','y','z']

    self.plot([center[0]],[center[1]],[center[2]],'ko')
    for v,i in enumerate(Perturb(np.zeros(3),length)):
      ends = i + center
      artist = Arrow3D(*zip(center,ends),**kwargs)
      self.text(ends[0],ends[1],ends[2],labels[v],color=[0.3,0.3,0.3])
      self.add_artist(artist)


  def equal(self):
    self.set_pbaspect([1.0,1.0,1.0])

    d = np.max([np.diff(self.get_xlim3d()),
                np.diff(self.get_ylim3d()),
                np.diff(self.get_zlim3d())])/2.0
    
    mx = np.mean(self.get_xlim3d())
    my = np.mean(self.get_ylim3d())
    mz = np.mean(self.get_zlim3d())
    
    self.auto_scale_xyz((mx-d,mx+d),(my-d,my+d),(mz-d,mz+d))


  def get_proj(self):
    """                                              
    Create the projection matrix from the current viewing position.                       
    elev stores the elevation angle in the z plane                                     
    azim stores the azimuth angle in the x,y plane                                             
    dist is the distance of the eye viewing point from the object                                                  
    point.                                       
    """
    relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

    xmin, xmax = self.get_xlim3d()/self.pbaspect[0]
    ymin, ymax = self.get_ylim3d()/self.pbaspect[1]
    zmin, zmax = self.get_zlim3d()/self.pbaspect[2]

    # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0  
    worldM = proj3d.world_transformation(xmin, xmax,
                                         ymin, ymax,
                                         zmin, zmax)

    # look into the middle of the new coordinates                                           
    R = np.array([0.5, 0.5, 0.5])

    xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
    yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
    zp = R[2] + np.sin(relev) * self.dist
    E = np.array((xp, yp, zp))

    self.eye = E
    self.vvec = R - E
    self.vvec = self.vvec / proj3d.mod(self.vvec)

    if abs(relev) > np.pi/2:
      # upside down                                                                                
      V = np.array((0, 0, -1))
    else:
      V = np.array((0, 0, 1))

    zfront, zback = -self.dist, self.dist

    viewM = proj3d.view_transformation(E, R, V)
    perspM = proj3d.persp_transformation(zfront, zback)
    M0 = np.dot(viewM, worldM)
    M = np.dot(perspM, M0)
    return M

    
if __name__ == '__main__':
  def f(points):
    return np.sqrt(np.sum(points**2,1))

  fig = plt.figure()
  ax = Axes3D(fig)
  xrot = np.pi/2.0001
  yrot = 0.0
  zrot = 0.0
  length = 2.0
  width = 2.0
  anchor = np.array([-1.0,0.0,1.0])

  ax.cross_section(f,anchor,
                   zrot,yrot,xrot,
                   length,width,cmap=tcm.slip_r,Nl=50,Nw=50,clim=(0.0,2.0),zorder=10)

  xrot = 0.01*np.pi
  yrot = 0.0
  zrot = 0.0
  length = 2.0
  width = 2.0
  anchor = np.array([-1.0,0.0,0.0])
  s = ax.cross_section(f,anchor,
                   zrot,yrot,xrot,
                   length,width,cmap=tcm.slip_r,Nl=50,Nw=50,zorder=0)
  plt.colorbar(s)
  ax.set_xlim(-3,3)
  ax.set_ylim(-3,3)
  ax.set_zlim(-3,3)


  plt.show()
  
