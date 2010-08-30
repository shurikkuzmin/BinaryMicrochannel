#!/usr/bin/python
import numpy
import matplotlib
matplotlib.use("GTKAgg")
import pylab
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def load_phase(dir):
   array=numpy.loadtxt("./"+dir+"/phase040000.dat")
   pylab.figure()
   pylab.imshow(array)

def load_file(dir):
    array=numpy.loadtxt("./"+dir+"/velocity040000.dat")
    
    ny=51
    nx=751
      
    #data=array[:, 2:3]
 
    #absnumpy=numpy.reshape(data,  (nz, ny))
    #pylab.imshow(absnumpy)
    pylab.figure()
    pylab.imshow(array)
    #pylab.savefig("test.eps", format="eps", dpi=200)
    
    return array
#    fig = pylab.figure()
#    ax = Axes3D(fig)
#    X = numpy.arange(0, nx)
#    Y = numpy.arange(0, ny)
#    X, Y = numpy.meshgrid(X, Y)
#    ax.plot_surface(X, Y,array, cmap=cm.jet)
    
def show_profile(array, coor):
    pylab.figure()
    pylab.plot(array[:, coor])

def compare_profiles(dir1, dir2, coor1, coor2):
    array_force=numpy.loadtxt("./"+dir1+"/velocity040000.dat")
    array_pressure=numpy.loadtxt("./"+dir2+"/velocity040000.dat")
    pylab.figure()
    pylab.plot(array_force[:, coor1], "g+")
    pylab.plot(array_pressure[:, coor2])
    
if __name__=="__main__":
    #print "Velocity in the center=",absnumpy[ny/2,nx/2]
    
    load_phase("Force/10/tmp2")
    load_phase("Pressure/10/tmp")

     #array=load_file("Pressure/14/tmp")
    #array=load_file("tmp2")
    #show_profile(array,550)
    #show_profile(array, 100)
    compare_profiles("Force/10/tmp2", "Pressure/10/tmp", 630, 550)
    pylab.show()
