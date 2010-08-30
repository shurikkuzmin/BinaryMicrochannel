#!/usr/bin/python
import os
from subprocess import call

if __name__=="__main__":
    print os.getcwd()
    #os.mkdir("temp")
    for i in range(1, 5):
        dir_temp=str(49*i+2)
        os.mkdir(dir_temp)
        os.chdir(dir_temp)
        os.mkdir("Pressure")
        os.mkdir("Force")
        subprocess.Popen('cp ../binary_simple.cpp Force/', shell=True)
        subprocess.Popen('cp ../binary_guo_pressure_bb_wall.cpp Pressure/', shell=True)

        os.chdir("..")
        
