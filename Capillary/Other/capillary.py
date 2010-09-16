#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    grid_ny=[314,210,158,127,106,98,85]
    grid_nx=[4681,3121,2341,1876,1561,1441,1246]
    thickness=[0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.15]
    capillary=[0.03, 0.05, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
    
    #os.mkdir("temp")
    for i in range(0, 7):
        dir_temp=str(grid_ny[i])
        os.mkdir(dir_temp)
        os.chdir(dir_temp)
        subprocess.Popen('cp ../binary_simple.cpp .', shell=True)
        os.mkdir("tmp2")
        cmd="sed -e 's/const int NY=[0-9]*;/const int NY="+str(grid_ny[i])+";/' -e 's/const int NX=[0-9]*;/const int NX="+str(grid_nx[i])+";/'"
        #print cmd
        sed=subprocess.Popen(cmd+" binary_simple.cpp > binary_simple_out.cpp", shell=True)
        sed.wait()
        print os.getcwd()
        comp=subprocess.Popen("g++ -O3 binary_simple_out.cpp -o main.out", shell=True)
        comp.wait()
        subprocess.Popen("./main.out "+str(thickness[i]*(grid_ny[i]-2))+" "+str(0.28812*capillary[i]/((grid_ny[i]-2)*(grid_ny[i]-2)))+" &",shell=True)
        os.chdir("..")

def Get_Zero(prof, thickness):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-0.5)/thickness

def Analyze_Simulations():
    print os.getcwd()
    grid_ny=[314,210,158,127,106,98,85]
    up=[127, 158, 210, 314]
    down=[85, 98, 106]
    lengths=[1700, 2000, 2400, 3000]
    
    for i in range(0, len(up)):
        dir_temp=str(up[i]) #str(49*i+2)
        os.chdir(dir_temp)
        #os.chdir("Force")
        pylab.figure()
        name="phase"+(1-(up[i]*1020)/100000)*"0"+str(1020*up[i])+".dat"
        array=numpy.loadtxt(name)
        pylab.imshow(array)
        #pylab.figure()
        #pylab.plot(array[:, lengths[i]])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(array[:, lengths[i]], up[i]-2)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("..")

def Analyze_Velocities():
    print os.getcwd()
    grid_ny=[314,210,158,127,106,98,85]
    up=[127, 158, 210, 314]
    down=[85, 98, 106]
    lengths=[1700, 2000, 2400, 3000]
    bulk=[750, 750, 500, 100]
    
    for i in range(0, len(up)):
        dir_temp=str(up[i]) #str(49*i+2)
        os.chdir(dir_temp)
        #print os.getcwd()
        #os.chdir("Force")
        name="velocity"+(1-(up[i]*1020)/100000)*"0"+str(1020*up[i])+".dat"
        array=numpy.loadtxt(name)
        prof=array[:,bulk[i]]
        pylab.figure()
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:, bulk[i]])
        print prof[len(prof)/2]
        os.chdir("..")


if __name__=="__main__":
    #Run_Simulations()
    #Analyze_Simulations()    
    Analyze_Velocities()
    pylab.show()
