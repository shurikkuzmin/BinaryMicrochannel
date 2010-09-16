#!/usr/bin/python
import os
import subprocess
import pylab
import numpy

if __name__=="__main__":
    dir=["Other","Results","SmallCapillary"]
    styles=["r+-.","b^-","go--"]
    for counter,name_dir in enumerate(dir):
        os.chdir(name_dir)
        file="capillary_number.dat"
        arr=numpy.loadtxt(file)
        
        capillary_aver=0.5*(arr[:,1]+arr[:,2])
        pylab.plot(capillary_aver,arr[:,0],styles[counter],markersize=14,linewidth=3)    
        os.chdir("..")
    pylab.xlabel("Capillary number",fontsize=16)
    pylab.ylabel("Thickness",fontsize=16)
    pylab.show()