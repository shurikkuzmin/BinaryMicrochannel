#!/usr/bin/python

import numpy
import os
import pylab


def Get_Zero(prof):   
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)

if __name__=="__main__":
    arr_sailfish=numpy.load("sailfish0003000.npz")
    phi_sailfish=arr_sailfish['phi']
    phi_my=numpy.loadtxt("phase003000.dat")
    prof_sailfish=phi_sailfish[:,400]
    prof_sailfish=prof_sailfish[1:len(prof_sailfish)-1]
    prof_my=phi_my[:,400]
    dist=len(prof_my)

    #bnd = np.argwhere(a['phi'] < 0)
    #midx = (np.min(bnd[:,1]) + np.max(bnd[:,1])) / 2
    #phi = a['phi'][:,midx]


    print prof_sailfish[2]," compare with ", prof_sailfish[dist-3]
    print prof_my[2], " compare with ", prof_my[dist-3]
    print Get_Zero(prof_sailfish)," compare with ", Get_Zero(prof_my)

    pylab.figure()
    pylab.imshow(arr_sailfish['phi'])
    pylab.figure()
    pylab.imshow(phi_my)



    pylab.figure()
    pylab.plot(prof_sailfish)
    pylab.plot(prof_my)
    pylab.show()
    
