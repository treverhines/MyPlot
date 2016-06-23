#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar
import matplotlib.colors

def pseudo_transparent_cmap(cmap,alpha):
  arr = cmap(np.arange(256))
  arr[:,[0,1,2]] = (1-alpha) + alpha*arr[:,[0,1,2]]
  new_cmap = matplotlib.colors.ListedColormap(arr)
  return new_cmap

def transparent_colorbar(*args,**kwargs):
  cbar = plt.colorbar(*args,**kwargs)
  if cbar.alpha is None:
    alpha = 1.0
  else:
    alpha = cbar.alpha
  cbar.set_alpha(1.0)

  cbar.cmap = pseudo_transparent_cmap(cbar.cmap,alpha)
  cbar.draw_all()
  cbar.solids.set_rasterized(True)
  return cbar
  
def overlay(top,bot):
  ''' 
  Parameters
  ----------
    rgba : (...,4) array

    base : (...,4) array or (...,3) array

  Returns
  -------
    out : (...,3) array 
  '''
  top = np.asarray(top)
  bot = np.asarray(bot)
  if bot.shape[-1] == 4:
    base_rgb = np.ones(bot.shape[:-1]+(3,))
    bot = overlay(bot,base_rgb)

  alpha = top[...,3]
  rgb = bot*(1-alpha) + alpha*top[...,[0,1,2]]
  return rgb

  

