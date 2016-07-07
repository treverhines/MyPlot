#!/usr/bin/env python
import numpy as np
import mayavi.mlab as mlab
import matplotlib.cm as cm
from colormaps.colormaps import viridis
import transform
import myplot.cm
from traits.api import HasTraits, Range, Instance, on_trait_change
from traitsui.api import View, Item, Group
from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor,MlabSceneModel

def cmap_to_rgba(cmap):
  x = np.linspace(0,1,256)
  rgba = np.array([cmap(i) for i in x])
  rgba *= 255
  rgba = np.array(rgba,dtype=int)
  return rgba

class VectorXSection:
  def __init__(self,f,
               f_args=None, 
               f_kwargs=None,
               Nl=20,Nw=20,
               base_square_x=(0,1),
               base_square_y=(0,1),
               transforms = (),
               scale_units=None):

    if f_args is None:
      f_args = ()

    if f_kwargs is None:
      f_kwargs = {}

    self._f = f
    self._f_args = f_args
    self._f_kwargs = f_kwargs
    self._Nl = Nl
    self._Nw = Nw
    self._base_square_x = base_square_x
    self._base_square_y = base_square_y
    self._transforms = transforms
    self._plots = None
    self.scale_units = scale_units  

  def set_f(self,f):
    self._f = f

  def set_f_kwargs(self,f_kwargs):
    self._f_kwargs = f_kwargs

  def set_f_args(self,f_args):
    self._f_args = f_args

  def set_base_square_x(self,a):
    self._base_square_x = a

  def set_base_square_y(self,a):
    self._base_square_y = a

  def add_transform(self,t):
    self._transforms += t,

  def draw(self,**kwargs):
    x = np.linspace(self._base_square_x[0] + 0.5*(self._base_square_x[1]-self._base_square_x[0])/self._Nl,
                    self._base_square_x[1] - 0.5*(self._base_square_x[1]-self._base_square_x[0])/self._Nl,
                    self._Nl)
    y = np.linspace(self._base_square_y[0] + 0.5*(self._base_square_y[1]-self._base_square_y[0])/self._Nw,
                    self._base_square_y[1] - 0.5*(self._base_square_y[1]-self._base_square_y[0])/self._Nw,
                    self._Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    for trans in self._transforms:
      p = np.concatenate((x[:,:,None],
                          y[:,:,None],
                          z[:,:,None]),
                          axis=-1)
      p = trans(p)
      pflat = np.reshape(p,(-1,3))
      x,y,z = pflat[:,0],pflat[:,1],pflat[:,2]
      value = self._f(pflat,*self._f_args,**self._f_kwargs)
      u,v,w = value[:,0],value[:,1],value[:,2]
      m = mlab.quiver3d(x,y,z,u,v,w,mode='arrow',color=(1.0,1.0,1.0),scale_factor=self.scale_units,
                        resolution=20,**kwargs)
      mlab.draw()
 

      if self._plots is None:
        self._plots = (m,trans),

      else:
        self._plots += (m,trans),

    return [i[0] for i in self._plots]

  def redraw(self):
    x = np.linspace(self._base_square_x[0],
                    self._base_square_x[1],self._Nl)
    y = np.linspace(self._base_square_y[0],
                    self._base_square_y[1],self._Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    for plot,trans in self._plots:
      p = np.concatenate((x[:,:,None],
                          y[:,:,None],
                          z[:,:,None]),
                          axis=-1)
      p = trans(p)
      pflat = np.reshape(p,(-1,3))
      x,y,z = pflat[:,0],pflat[:,1],pflat[:,2]
      value = self._f(pflat,*self._f_args,**self._f_kwargs)
      u,v,w = value[:,0],value[:,1],value[:,2]
      plot.mlab_source.set(x=x,y=y,z=z,
                           u=u,v=v,w=w)

      if self._plots is None:
        self._plots += (m,trans),

    return [i[0] for i in self._plots]

  def view(self):
    mlab.show()


class XSection:
  def __init__(self,f,
               f_args=None, 
               f_kwargs=None,
               Nl=400,Nw=400,
               base_square_x=(0,1),
               base_square_y=(0,1),
               transforms = (),
               clim=None,
               cmap=viridis):

    if f_args is None:
      f_args = ()

    if f_kwargs is None:
      f_kwargs = {}

    self._f = f
    self._f_args = f_args
    self._f_kwargs = f_kwargs
    self._Nl = Nl
    self._Nw = Nw
    self._clim = clim
    self._rgba = cmap_to_rgba(cmap)
    self._base_square_x = base_square_x
    self._base_square_y = base_square_y
    self._transforms = transforms
    self._plots = None

  def set_f(self,f):
    self._f = f

  def set_f_kwargs(self,f_kwargs):
    self._f_kwargs = f_kwargs

  def set_f_args(self,f_args):
    self._f_args = f_args

  def set_base_square_x(self,a):
    self._base_square_x = a

  def set_base_square_y(self,a):
    self._base_square_y = a

  def add_transform(self,t):
    self._transforms += t,

  def draw(self,**kwargs):
    x = np.linspace(self._base_square_x[0],
                    self._base_square_x[1],self._Nl)
    y = np.linspace(self._base_square_y[0],
                    self._base_square_y[1],self._Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    for trans in self._transforms:
      p = np.concatenate((x[:,:,None],
                          y[:,:,None],
                          z[:,:,None]),
                          axis=-1)
      p = trans(p)
      pflat = np.reshape(p,(-1,3))
      c = self._f(pflat,*self._f_args,**self._f_kwargs)
      c = np.reshape(c,np.shape(p)[:-1])
      if self._clim is None:
        self._clim = (np.min(c),np.max(c))

      m = mlab.mesh(p[:,:,0],p[:,:,1],p[:,:,2],
                    scalars=c,
                    vmin=self._clim[0],
                    vmax=self._clim[1],**kwargs)

      m.module_manager.scalar_lut_manager.lut.table = self._rgba
      mlab.colorbar()
      mlab.draw()

      if self._plots is None:
        self._plots = (m,trans),

      else:
        self._plots += (m,trans),

    return [i[0] for i in self._plots]

  def redraw(self):
    x = np.linspace(self._base_square_x[0],
                    self._base_square_x[1],self._Nl)
    y = np.linspace(self._base_square_y[0],
                    self._base_square_y[1],self._Nw)
    x,y = np.meshgrid(x,y)
    z = 0.0*x
    for plot,trans in self._plots:
      p = np.concatenate((x[:,:,None],
                          y[:,:,None],
                          z[:,:,None]),
                          axis=-1)
      p = trans(p)
      pflat = np.reshape(p,(-1,3))
      c = self._f(pflat,*self._f_args,**self._f_kwargs)
      c = np.reshape(c,np.shape(p)[:-1])
      plot.mlab_source.set(x=p[:,:,0],
                           y=p[:,:,1],
                           z=p[:,:,2],
                           scalars=c)
      if self._plots is None:
        self._plots += (m,trans),

    return [i[0] for i in self._plots]

  def view(self):
    mlab.show()


if __name__ == '__main__':
  import numpy as np
  from numpy import arange, pi, cos, sin

  from traits.api import HasTraits, Range, Instance, \
    on_trait_change, Int
  from traitsui.api import View, Item, Group

  from mayavi.core.api import PipelineBase
  from mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel

  import transform

  x = np.random.random((10))
  y = np.random.random((10))
  z = np.random.random((10))

  def foo(x):
    return x


  u = np.random.random((10))
  v = np.random.random((10))
  w = np.random.random((10))
  mlab.quiver3d(x,y,z,u,v,w,mode='arrow')
  mlab.show()

  dphi = pi/1000.
  phi = arange(0.0, 2*pi + 0.5*dphi, dphi, 'd')

  def f(p,farg1,farg2):
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    t = np.sqrt(farg1*x**2 + farg2*y**2)
    return t

  t = transform.identity()
  print('start1')
  m = VectorXSection(foo,transforms=[t])
  print('end1')
  m.draw()
  m.view()

  class InteractiveXSection(HasTraits):
    # define range of f_arg values             
    farg1 = Range(0, 30, .6)
    farg2 = Range(0, 30, .11)
    scene = Instance(MlabSceneModel, ())
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                Group('farg1', 'farg2'),
                resizable=True)

    def __init__(self):
      # initiate cross section instance            
      self.xsection = XSection(f,
                               Nl = 20,
                               Nw = 20,clim=(0,0.1),
                               f_args=(self.farg1,self.farg2))
      # add cross sections                                              
      T = transform.identity()
      T += transform.point_stretch([2.0,2.0,1.0])
      T += transform.point_translation([-1.0,0.0,0.0])
      self.xsection.add_transform(T)
      T += transform.point_rotation_x(np.pi/2.0)
      self.xsection.add_transform(T)
      T += transform.point_rotation_x(np.pi/2.0)
      self.xsection.add_transform(T)
      T += transform.point_rotation_x(np.pi/2.0)
      self.xsection.add_transform(T)
      # run HasTraits initiation routine                                                               
      HasTraits.__init__(self)
      # When the scene is activated, or when the parameters are changed, we   
      # update the plot.                            

    @on_trait_change('farg1,farg2,scene.activated')
    def update_plot(self):
      self.xsection.set_f_args((self.farg1,self.farg2))
      if self.xsection._plots is None:
        self.xsection.draw()

      else:
        self.xsection.redraw()
    
  my_model = InteractiveXSection()
  my_model.configure_traits()
  

