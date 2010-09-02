#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    #os.mkdir("temp")
    for i in range(1, 5):
        dir_temp=str(49*i+2)
        os.mkdir(dir_temp)
        os.chdir(dir_temp)
        os.mkdir("Pressure")
        os.mkdir("Force")
        subprocess.Popen('cp ../../binary_simple.cpp Force/', shell=True)
        subprocess.Popen('cp ../../binary_guo_pressure_bb_walls.cpp Pressure/', shell=True)
        os.chdir("Force")
        os.mkdir("tmp2")
        cmd="sed -e 's/const int NY=[0-9]*;/const int NY="+str(49*i+2)+";/' -e 's/const int NX=[0-9]*;/const int NX="+str(15*49*i+1)+";/'"
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

def Analyze_Simulations():
    print os.getcwd()
    from scipy.interpolate import UnivariateSpline
    
    for i in range(1, 4):
        dir_temp=str(49*i+2)
        os.chdir(dir_temp)
        os.chdir("Force")
        pylab.figure()
        name="phase"+(1-(i*40000)/100000)*"0"+str(40000*i)+".dat"
        array=numpy.loadtxt(name)
        pylab.imshow(array)
        pylab.figure()
        pylab.plot(array[:, 600*i])
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("../..")
if __name__=="__main__":
    Analyze_Simulations()    
    pylab.show()
