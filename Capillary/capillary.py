#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    grid_ny=[314,210,158,127,106,98,85]
    grid_nx=[4681,3121,2341,1876,1561,1441,1246]
    #os.mkdir("temp")
    for i in range(1, 5):
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
        comp=subprocess.Popen("g++ -O3 binary_simple_out.cpp -o main.out")
        comp.wait()
        subprocess.Popen("./main.out "+str(6*i)+" "+str(i)+" &")
        os.chdir("..")
        os.chdir("Pressure")
        os.mkdir("tmp")
        sed=subprocess.Popen(cmd+" binary_guo_pressure_bb_walls.cpp > binary_guo_pressure_bb_walls_out.cpp",  shell=True)
        sed.wait()
        comp=subprocess.Popen("g++ -O3 binary_guo_pressure_bb_walls_out.cpp -o main.out")
        comp.wait()
        subprocess.Popen("./main.out "+str(6*i)+" "+str(i)+" &")
        os.chdir("..")
        os.chdir("..")

def Get_Zero(prof, factor):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-0.5)/(49*factor)

def Analyze_Simulations():
    print os.getcwd()
    
    for i in range(1, 5):
        dir_temp=str(49*i+2)
        os.chdir(dir_temp)
        os.chdir("Force")
        pylab.figure()
        name="phase"+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:, 600*i])
        pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(array[:, 600*i], i)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("../..")

def Analyze_Velocities():
    print os.getcwd()
    
    for i in range(1, 5):
        dir_temp=str(49*i+2)
        os.chdir(dir_temp)
        os.chdir("Force")
        name="velocity"+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        prof=array[:,100]
        pylab.figure()
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:, 100])
        print prof[len(prof)/2]
        os.chdir("../..")


if __name__=="__main__":
    Analyze_Simulations()    
    #Analyze_Velocities()
    pylab.show()
