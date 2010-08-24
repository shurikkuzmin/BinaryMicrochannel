#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if __name__=="__main__":
    array=numpy.loadtxt("./tmp/phase000300.dat")
    
    ny=31
    nx=201
      
    #data=array[:, 2:3]
 
    #absnumpy=numpy.reshape(data,  (nz, ny))
    #pylab.imshow(absnumpy)
    pylab.imshow(array)
    fig = pylab.figure()
    ax = Axes3D(fig)
    X = numpy.arange(0, nx)
    Y = numpy.arange(0, ny)
    X, Y = numpy.meshgrid(X, Y)
    ax.plot_surface(X, Y,array, cmap=cm.jet)
    #print "Velocity in the center=",absnumpy[ny/2,nx/2]

    pylab.show()
