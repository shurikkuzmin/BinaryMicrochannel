#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    ny=[102,127,152,177,202,227]
    nx=[1501,1876,2251,2626,3001,3376]
    force_init=6e-6/4;
    init_width=6
    #os.mkdir("temp")
    for i in range(0, 1):
        dir_temp=str(ny[i])
        os.mkdir(dir_temp)
        subprocess.call(['cp','bmc.py',dir_temp+"/"])
        os.chdir(dir_temp)
        ratio=float(ny[i]-2)/100.0
        width=int(init_width*ratio)
        subprocess.call(['./bmc.py','--lat_nx='+str(nx[i]),'--lat_ny='+str(ny[i]),'--bc_wall=halfbb',\
            '--force='+str(force_init/(ratio*ratio)),'--batch','--every=50000','--max_iters='+str(int(200000*ratio)+1),'--output=grid','--iwidth='+str(width),'--bc_wall_grad_order=1'])
        
        os.chdir("..")

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    return (zero-0.5)/(len(prof)-2)

    
def Analyze_Simulations():
    print os.getcwd()
        
    for i in range(0,400000):
        dir_temp="Proper Grid/"+str(49*i+2)
        os.chdir(dir_temp)
        #os.chdir("Force")
        pylab.figure()
        
        
def Analyze_Simulations():
    print os.getcwd()
    
    
    lengths=[]
    thicknesses=[]
    velocities=[]
    reynolds=[]
    for i in range(0,800000,20000):
        #os.chdir("Force")
        name="steady"+(6-len(str(i)))*"0"+str(i)+".npz"
        print name        
        array=numpy.load(name)
        phase=array['phi']
        vel=array['v'][0]
        dims=phase.shape
        center=phase[dims[0]/2,:]
           
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[1]
            z1=numpy.max(numpy.where(center>0.0))
        print z1,z2
        lengths.append(z2-z1)
        thicknesses.append(Get_Zero(phase[:,((z1+z2)/2)%dims[1]]))
        velocities.append(vel[dims[0]/2,z2%dims[1]])        
        reynolds.append(vel[dims[0]/2,z2%dims[1]]*200.0/(2.0/3.0))        
        #pylab.figure()
        #pylab.imshow(array['phi'])
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(array[:, 520*i], i)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)
    print "Lengths=",lengths
    print "Thicknesses=",thicknesses
    print "Velocities=",velocities
    print "Reynolds=",reynolds

if __name__=="__main__":
    Analyze_Simulations()    
    #Analyze_Velocities()
    #Run_Simulations()
    #pylab.show()
