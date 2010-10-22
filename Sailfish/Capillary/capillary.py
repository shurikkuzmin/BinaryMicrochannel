#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    capillary=[0.03,0.05,0.08,1.0,0.2,0.4,0.6,0.8,1.0]
    capillary_str=["3","5","8","10","20","40","60","80","100"]
    width=[0.04,0.06,0.08,0.1,0.12,0.13,0.15]
    force_init=6e-6/16;
    init_width=12
    #os.mkdir("temp")
    for i in range(0, len(capillary)):
        dir_temp=capillary_str[i]
        os.mkdir(dir_temp)
        subprocess.call(['cp','bmc.py',dir_temp+"/"])
        os.chdir(dir_temp)
        ratio=capillary[i]/0.05
        force=force_init*ratio
        width_value=int(width[i]*200)
        subprocess.call(['./bmc.py','--bc_wall=halfbb','--force='+str(force),'--batch','--every=50000','--max_iters=200000','--output=capillary','--iwidth='+str(width_value)])
        
        os.chdir("..")

def Get_Zero(prof):
    zero=0
    #pylab.figure()
    #pylab.plot(prof)
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-0.5)/(len(prof)-2)

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
    #Analyze_Velocities()
    Run_Simulations()
    #pylab.show()
