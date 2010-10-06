#!/usr/bin/python
import os
import subprocess
import pylab
import numpy

def Get_Zero(prof):
    zero=0
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)


if __name__=="__main__":
    dir=["Other","Results","SmallCapillary"]
    styles=["r+-.","b^-","go--"]
    
    for counter in range(-10, 12, 2):
        name_dir='grad'+str(counter)
        os.chdir(name_dir)
        file="phase200000.dat"
        arr=numpy.loadtxt(file)
        #pylab.plot(capillary_aver,arr[:,0],styles[counter],markersize=14,linewidth=3)    
        pylab.figure()
        #pylab.imshow(arr)
        prof=arr[:, 2600]
        print Get_Zero(prof)
        pylab.plot(prof)
        os.chdir("..")
    
    #pylab.xlabel("Capillary number",fontsize=16)
    #pylab.ylabel("Thickness",fontsize=16)
    pylab.show()
