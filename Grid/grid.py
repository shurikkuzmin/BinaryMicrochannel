#!/usr/bin/python
import os
import subprocess

if __name__=="__main__":
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
        
