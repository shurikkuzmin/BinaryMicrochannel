#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
import math
 
def Run_Simulations():
    print os.getcwd()
    capillary=[0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8,1.0]
    capillary_str=["3","5","8","10","20","40","60","80","100"]
    width=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16,0.17]
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
        subprocess.call(['./bmc.py','--bc_wall=halfbb','--force='+str(force),'--batch','--every=50000','--max_iters=200001','--output=capillary','--iwidth='+str(width_value)])
        
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
    from numpy import genfromtxt
    print os.getcwd()
    capillary_theor=[0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8,1.0]
    capillary_str=["3","5","8","10","20","40","60","80","100"]
    width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16,0.17]

    #exam=[2000,2200,2600,2500,1100,950,600,350,350] Old values
    exam=[2000,2100,2200,2400,100,1600,200,1600,350]
    good=[0,1,2,3,4,5,6,7,8]
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    
    #read the giavedoni data (not precise though)
    #giavedoni=genfromtxt("planarcasesolution.csv",delimiter=',',dtype=float)[1:]
    
    #read old data
    capillary_not_viscous=numpy.loadtxt("../Capillary/capillary.dat")


    for i in range(0, len(capillary_theor)):
        dir_temp="Results/"+capillary_str[i]
        #ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="capillary200000.npz"
        array=numpy.load(name)
        prof=array['phi'][:,exam[i]]
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        #pylab.imshow(array['phi'])
    
        if i in good:
            #pylab.plot(prof)
            widths.append(Get_Zero(prof))
            vel=array['v'][0]
            vel_prof=vel[:,exam[i]]
            velocities.append(vel_prof[len(vel_prof)/2])
        
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)

        os.chdir("../..")
    
    fig=pylab.figure()
    capillaries=numpy.array(velocities)*(4.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    print "Widths=",widths
    print "Capillaries=",capillaries
    print "Velocities",velocities
    
    pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)
    pylab.loglog(capillary_not_viscous[:,0],capillary_not_viscous[:,1],"bD-",linewidth=3,markersize=10)
    #pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    #pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    pylab.xlim(0.02,2.5)
    pylab.ylim(ymin=0.01)
    numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend([r'''$\frac{\mu_{liq}}{\mu_{gas}}=20$''',r'''$\frac{\mu_{liq}}{\mu_{gas}}=10$'''])
    #pylab.xlim(xmax=15)
    pylab.savefig("capillaries_viscous.eps",format="EPS",dpi=300)


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
    Analyze_Simulations()    
    #Analyze_Velocities()
    #Run_Simulations()
    pylab.show()
