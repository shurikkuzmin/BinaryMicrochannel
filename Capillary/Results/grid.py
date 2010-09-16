#!/usr/bin/python
import os
import subprocess
import pylab
import numpy

def Get_Zero(prof, thickness):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-0.5)/thickness

def Analyze_Simulations():
    up=[1,2]
    down=[3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    print os.getcwd()
    analyze=[2600,350]
    for i in range(0, len(up)):
        dir_temp=str(up[i])
        os.chdir(dir_temp)
        #os.chdir("Force")
        pylab.figure()
        name="phase200000.dat" #+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        pylab.imshow(array)
        pylab.figure()
        prof=array[:,analyze[i]]
        pylab.plot(prof)
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(prof, len(prof))
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("..")

def Analyze_Velocities():
    up=[1,2]
    down=[3,4,5,6,7,8,9,10,11,12,13,14,15,16]

    analyze=[2600,350]
    bulk=[1250,500]
    print os.getcwd()
    
    for i in range(0,len(up)):
        dir_temp=str(up[i])
        os.chdir(dir_temp)
        #os.chdir("Force")
        name="velocity200000.dat" #"velocity"+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        #prof=array[:,analyze[i]]
        prof=array[:,bulk[i]]
        pylab.figure()
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:,analyze[i]])
        print prof[len(prof)/2]
        os.chdir("..")


if __name__=="__main__":
    #Analyze_Simulations()    
    Analyze_Velocities()
    pylab.show()
