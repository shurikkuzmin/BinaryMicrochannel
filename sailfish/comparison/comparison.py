#!/usr/bin/python

import numpy
import os
import pylab

if __name__=="__main__":
    arr_sailfish=numpy.load("sailfish0003000.npz")
    phi_sailfish=arr_sailfish['phi']
    phi_my=numpy.loadtxt("phase003000.dat")
    prof_sailfish=phi_sailfish[:,400]
    prof_sailfish=prof_sailfish[1:len(prof_sailfish)-1]
    prof_my=phi_my[:,400]
    dist=len(prof_my)

    print prof_sailfish[2]," compare with ", prof_sailfish[dist-3]
    print prof_my[2], " compare with ", prof_my[dist-3]

    pylab.figure()
    pylab.imshow(arr_sailfish['phi'])
    pylab.figure()
    pylab.imshow(phi_my)



    pylab.figure()
    pylab.plot(prof_sailfish)
    pylab.plot(prof_my)
    pylab.show()
    
