#!/usr/bin/python
import os
import subprocess
import pylab
import numpy
 
def Run_Simulations():
    print os.getcwd()
    gradient=numpy.arange(-1.0,1.2,0.2)
    str_gradient=['-10','-8','-6','-4','-2','0','2','4','6','8','10']
    for i in range(0, len(gradient)):
        dir_temp="Grad"+str_gradient[i]
        os.mkdir(dir_temp)
        subprocess.call(['cp','binary_microchannel.py',dir_temp+"/"])
        os.chdir(dir_temp)
        subprocess.call(['./binary_microchannel.py','--bc_wall=halfbb','--phase_grad='+str(gradient[i]),'--every=20000','--max_iters=200000','--output=grad','--batch'])
        os.chdir("..")
  
def Get_Zero(prof):
    zero=0
    for counter in range(0, len(prof)/2):
        if prof[counter]>=0 and prof[counter+1]<0:
            zero=-(prof[counter]*(counter+1)-prof[counter+1]*counter)/(prof[counter+1]-prof[counter])
    print (zero-1.5)/(len(prof)-4)

def Analyze_Simulations():
    print os.getcwd()
    gradient=[-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
    str_gradient=['-10','-8','-6','-4','-2','0','2','4','6','8','10']
    style=["bo","rH","g--","k-.","c<-.","y>","bs","g^","rD"]
    fig=pylab.figure()
    for i in range(0, len(str_gradient)-2):
        dir_temp="Results/Grad"+str_gradient[i]
        os.chdir(dir_temp)
        name="grad180000.npz"
        array=numpy.load(name)
        #fig=pylab.figure(figsize=(11,1.5))
        #pylab.imshow(array['phi'])
        pylab.plot(array['phi'][0:22,2000],style[i],linewidth=3,markersize=10)
        #pylab.yticks([0,25,50,75,100])
        #pylab.savefig("../../initfinish"+str(width[i])+".eps",format="EPS",dpi=70)
        
        #fig_init=pylab.figure(figsize=(11,1.5))
        #array_init=numpy.load("init000000.npz")
        #pylab.imshow(array_init['phi'])
        #pylab.yticks([0,25,50,75,100])
        #pylab.savefig("../../initbegin"+str(width[i])+".eps",format="EPS",dpi=70)
        
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        Get_Zero(array['phi'][:, 2000])
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)
        
        os.chdir("../..")
    fig.subplots_adjust(left=0.15)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$\phi$''',fontsize=30)
    pylab.legend([r'''$\partial_n\phi=$'''+r'''$'''+str(x)+r'''$''' for x in gradient[:-2]])
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.savefig("phase_grad_profiles.eps",format="EPS",dpi=300)

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
