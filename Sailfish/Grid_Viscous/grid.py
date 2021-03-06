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
    for i in range(0, len(ny)):
        dir_temp=str(ny[i])
        os.mkdir(dir_temp)
        subprocess.call(['cp','bmc.py',dir_temp+"/"])
        os.chdir(dir_temp)
        ratio=float(ny[i]-2)/100.0
        width=int(init_width*ratio)
        subprocess.call(['./bmc.py','--lat_nx='+str(nx[i]),'--lat_ny='+str(ny[i]),'--bc_wall=halfbb',\
            '--force='+str(force_init/(ratio*ratio)),'--batch','--every=50000','--max_iters='+str(int(200000*ratio)+1),'--output=grid','--iwidth='+str(width)])
        
        os.chdir("..")

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-1.5)/(len(prof)-4)

def Analyze_Simulations():
    print os.getcwd()
    ny=[102,127,152,177,202,227]
    nx=[1501,1876,2251,2626,3001,3376]
    style=["bo","rH","c<","y>","bs","g^"]
    color=["b","r","c","y","b","g"]
    style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    fig=pylab.figure()


    for i in range(0, len(ny)):
        dir_temp="Results/"+str(ny[i])
        ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="grid"+str(200000+i*50000)+".npz"
        array=numpy.load(name)
        prof=array['phi'][:,int(1400*ratio)]
        x=numpy.arange(0.0,float(ny[i]))/ratio
        #pylab.imshow(array['phi'])
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("../..")
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$\phi$''',fontsize=30)
    labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend(labels)
    pylab.xlim(xmax=15)
    pylab.savefig("norm_grid_profs.eps",format="EPS",dpi=300)


def Analyze_Velocities():
    print os.getcwd()
    ny=[102,127,152,177,202,227]
    nx=[1501,1876,2251,2626,3001,3376]

    for i in range(0, len(ny)):
        dir_temp="Results/"+str(ny[i])
        ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        name="grid"+str(200000+50000*i)+".npz"
        array=numpy.load(name)
        ux=array['v'][0]
        #print ux.shape
        prof=ux[:,int(1400*ratio)]
        #pylab.figure()
        #pylab.plot(prof)
        print prof[len(prof)/2]
        os.chdir("../..")


if __name__=="__main__":
    #Analyze_Simulations()    
    Analyze_Velocities()
    #Run_Simulations()
    #pylab.show()
