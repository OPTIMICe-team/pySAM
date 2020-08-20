#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:17:59 2020

@author: dori
"""

#%gui qt # activate this magic if using ipython, or set --gui=qt at invocation

import numpy as np
from matplotlib import colors
from mayavi import mlab

data = np.loadtxt('../67_0.992257_1.96799e-05.dat')#, dtype=int)

def visualize(X, bgcolor=(1,1,1), fgcolor=(.8,.8,.8)):
	"""Visualize the aggregate using Mayavi.
	Args:
		bgcolor: Background color for the Mayavi scene.
		fgcolor: Foreground color for the Mayavi scene.
	"""

	color_list = [colors.colorConverter.to_rgb(c) for c in [
		"#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
		"#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
		"#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
	]]

	# local import as this can take a while
	#from mayavi import mlab
	
	mlab.figure(bgcolor=bgcolor, fgcolor=fgcolor)
	ident = np.unique(X[:, 3])
	for i, ID in enumerate(ident):
		x = X[X[:, 3]==ID]
		print(x)
		mlab.points3d(x[:,0], x[:,1], x[:,2],
				      color=color_list[i%len(color_list)],
				      mode="cube", scale_factor=1)
	mlab.savefig('figure.png')

visualize(data)
