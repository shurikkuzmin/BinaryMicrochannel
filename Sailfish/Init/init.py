#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    width=[12,14,16,18,20]
    force_init=6e-6/4;
    #os.mkdir("temp")
    for i in range(0, len(width)):
        dir_temp=str(width[i])
        os.mkdir(dir_temp)
        subprocess.call(['cp','binary_microchannel.py',dir_temp+"/"])
        os.chdir(dir_temp)
        subprocess.call(['./binary_microchannel.py','--bc_wall=halfbb','--init_width='+dir_temp,'--every=20000','--max_iters=200000','--output=init'])
        os.chdir("..")
  
def Get_Zero(prof):
    zero=0
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-0.5)/(len(prof)-2)

def Analyze_Simulations():
    print os.getcwd()
    width=[12,14,16,18,20]
    for i in range(0, len(width)):
        dir_temp="Results/"+str(width[i])
        os.chdir(dir_temp)
        name="init180000.npz"
        array=numpy.load(name)
        fig=pylab.figure(figsize=(11,1.5))
        pylab.imshow(array['phi'],cmap="gray",extent=[0.0,15.0,0.0,1.0])
        pylab.yticks([0.0,0.5,1.0])
        #pylab.axis([0.0,15.0,0.0,1.0])
        pylab.savefig("../../initfinish"+str(width[i])+".eps",format="EPS",dpi=70)
        
        fig_init=pylab.figure(figsize=(11,1.5))
        array_init=numpy.load("init000000.npz")
        pylab.imshow(array_init['phi'],cmap="gray",extent=[0.0,15.0,0.0,1.0])
        #pylab.axis([0.0,15.0,0.0,1.0])
        #pylab.xlim(0.0,15.0)
        #pylab.ylim(0.0,1.0)

        pylab.yticks([0.0,0.5,1.0])
        pylab.savefig("../../initbegin"+str(width[i])+".eps",format="EPS",dpi=70)
        
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(array['phi'][:, 1300])
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)
        os.chdir("../..")
        

def Analyze_Velocities():
    print os.getcwd()
    
    for i in range(1, 5):
        dir_temp="Proper Grid/"+str(49*i+2)
        os.chdir(dir_temp)
        #os.chdir("Force")
        name="velocity"+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        prof=array[:,520*i]
        pylab.figure()
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:, 520*i])
        print prof[len(prof)/2]
        os.chdir("../..")


if __name__=="__main__":
    Analyze_Simulations()    
    #Analyze_Velocities()
    #Run_Simulations()
    pylab.show()
