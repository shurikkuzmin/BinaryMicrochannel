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
    capillary_theor=[0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8]
    capillary_str=["3","5","8","10","20","40","60","80"]
    width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16]

    exam=[2000,2200,2600,2500,1100,950,600,350]
    good=[0,1,2,3,4,5,6,7]
    #style=["bo","rH","c<","y>","bs","g^"]
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    widths=[]
    velocities=[]
    velocities_real=[]
    re_real=[]
    widths_real=[]
    #read the giavedoni data (not precise though)
    giavedoni=genfromtxt("planarcasesolution.csv",delimiter=',',dtype=float)[1:]


    for i in range(0, len(capillary_theor)):
        dir_temp="Results/"+capillary_str[i]
        #ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="capillary200000.npz"
        array=numpy.load(name)
        prof=array['phi'][:,exam[i]]
        dims=array['phi'].shape
        
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        #pylab.imshow(array['phi'])
    
        if i in good:
            #pylab.plot(prof)
            widths.append(Get_Zero(prof))
            vel=array['v'][0]
            center=array['phi'][dims[0]/2,:]
           
            z1 = numpy.min(numpy.where(center < 0.0))
            z2 = numpy.max(numpy.where(center < 0.0))
            if z1==0:
                z2=numpy.min(numpy.where(center>0.0))+dims[1]
                z1=numpy.max(numpy.where(center>0.0))
            print z1,z2
            
            prof_real=array['phi'][:,((z1+z2)/2)%dims[1]]
            widths_real.append(Get_Zero(prof_real))     

            vel_prof=vel[:,exam[i]]
            velocities.append(vel_prof[len(vel_prof)/2])
            velocities_real.append(vel[dims[0]/2,z2%dims[1]])
            re_real.append(vel[dims[0]/2,z2%dims[1]]*dims[0]/(2.0/3.0))
        
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)

        os.chdir("../..")
    
    fig=pylab.figure()
    capillaries=numpy.array(velocities)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)
    capillaries_real=numpy.array(velocities_real)*(2.0/3.0)/math.sqrt(8.0*0.04*0.04/9.0)

    print "Widths=",widths
    print "Widths_real=",widths_real
    print "Capillaries=",capillaries
    print "Capillaries_real=",capillaries_real    
    print "Velocities=",velocities
    print "Velocities_real=",velocities_real
    print "Reynolds=",re_real

    
    pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    pylab.loglog(capillaries,widths,"go-",linewidth=3,markersize=10)

    pylab.xlim(0.02,1.5)
    pylab.ylim(ymin=0.01)
    numpy.savetxt("capillary.dat",zip(capillaries_real,widths_real))
    
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend(["Giavedoni","Heil","This work"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("capillaries_comparison.eps",format="EPS",dpi=300)

    #Another figure
    fig=pylab.figure()    
    pylab.loglog(giavedoni[:,0],giavedoni[:,1]/2.0,"bD-",linewidth=3,markersize=10)
    pylab.loglog(capillary_theor,width_theor,"ys--",linewidth=3,markersize=10)
    pylab.loglog(capillaries_real,widths,"go-",linewidth=3,markersize=10)

    pylab.xlim(0.02,1.5)
    pylab.ylim(ymin=0.01)
    numpy.savetxt("capillary.dat",zip(capillaries,widths))
    
    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$Ca$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    pylab.legend(["Giavedoni","Heil","This work"],loc=4)
    pylab.savefig("capillaries_comparison_real.eps",format="EPS",dpi=300)

    
def Analyze_Bubble():
    from numpy import genfromtxt
    print os.getcwd()
    capillary_theor=numpy.array([0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8])
    capillary_str=numpy.array(["3","5","8","10","20","40","60","80"])
    width_theor=[0.04,0.06,0.08,0.1,0.12,0.13,0.15,0.16]

    exam=[2000,2200,2600,2500,1100,950,600,350]
    good=[0,1,2,3,4,5,6,7]
    good_value=[1,2,4]
    style=["b-.","r--","c."]
    #labels=[r'''$Ca='''+ca+r'''$''' for ca in numpy.array(capillary_theor[good_value],dtype=str)]
    labels=[]    
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    low=[1509,1865,381]
    high=[2554,2900,1427]
    
    #ux=array['v'][0]
    #center=array['phi'][ny[i]/2,:]

    
    
    widths=[]
    velocities=[]
    fig=pylab.figure()
    for counter,value in enumerate(good_value):
        dir_temp="Results/"+capillary_str[value]
        #ratio=float(ny[i]-2)/100.0
        os.chdir(dir_temp)
        #os.chdir("Force")
        #pylab.figure()
        name="capillary200000.npz"
        array=numpy.load(name)
        thicknesses=[]
        dims=array['phi'].shape        
        ux=array['v'][0]
        center=array['phi'][dims[0]/2,:]
           
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+nx[i]
            z1=numpy.max(numpy.where(center>0.0))
        print z1,z2

        for coorz in range(z1,z2):
        #for coor in range(low[counter],high[counter]):
            coor=coorz%dims[1]
            prof=array['phi'][:,coor]
            thicknesses.append(Get_Zero(prof))
        #prof=array['phi'][:,exam[i]]
        #x=numpy.arange(0.0,float(ny[i]))/ratio
        
        #pylab.imshow(array['phi'])
    
        #if i in good:
            #pylab.plot(prof)
        #    widths.append(Get_Zero(prof))
        #    vel=array['v'][0]
        #    vel_prof=vel[:,exam[i]]
        #    velocities.append(vel_prof[len(vel_prof)/2])
        
        #pylab.plot(x[0:20],prof[0:20],color[i],linewidth=3)
        #pylab.plot(x,prof,style_diff[i],markersize=10,linewidth=3)
        #pylab.figure()
        #pylab.plot(array[:, 520*i])
        #pylab.savefig("grid_phase_prof_"+str(49*i)+".eps", dpi=300)
        #Get_Zero(prof)
        #extrapolator=UnivariateSpline(array[0:(49*i+2)/2, 600*i], numpy.arange(0, (49*i+2)/2),  k=2)
        #print extrapolator(0)
        
        #pylab.figure()
                
        print thicknesses[((z1+z2)/2)%dims[1]-z1]        
        print numpy.std(thicknesses[300:-300])/thicknesses[((z1+z2)/2)%dims[1]-z1]        
        delta_x=15.0/(len(thicknesses)-1)
        x=delta_x*numpy.arange(0,len(thicknesses))*float(len(thicknesses))/3000.0
        pylab.plot(x,thicknesses,style[counter],linewidth=3)
        labels.append(r'''$Ca='''+str(ux[dims[0]/2,z2%dims[1]]*2.0/3.0/math.sqrt(8.0*0.04*0.04/9.0))[0:4]+r'''$''')
        os.chdir("../..")
    
    pylab.ylim(ymin=0.0,ymax=0.5)
    pylab.xlim(xmax=5.1)

    leg=pylab.legend(labels)
    legtext = leg.get_texts() # all the text.Text instance in the legend
    for text in legtext:
        text.set_fontsize(20) 
    #set(ltext, fontsize='large') # the legend text fontsize

    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(["Simulations","Giavedoni","Heil"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("bubble_length.eps",format="EPS",dpi=300)

def smoothListGaussian(list,strippedXs=False,degree=2):  
    window=degree*2-1  
    weight=numpy.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(numpy.exp((4*(frac))**2))  
        weightGauss.append(gauss)  
    weight=numpy.array(weightGauss)*weight  
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
    return smoothed  

def smoothList(list,strippedXs=False,degree=5):  
    smoothed=[0]*(len(list)-degree+1)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(list[i:i+degree])/float(degree)  
    return smoothed  

def Analyze_Sehgal_Bubble():
    from numpy import genfromtxt
    print os.getcwd()
    style=["b-.","r--","c."]
    #labels=[r'''$Ca='''+ca+r'''$''' for ca in numpy.array(capillary_theor[good_value],dtype=str)]
    labels=[r'''$Ca='''+str(0.0113)+r'''$''',r'''$Ca='''+str(0.07)+r'''$''']   
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #fig=pylab.figure()
    #low=[1509,1865,381]
    #high=[2554,2900,1427]
    
    #ux=array['v'][0]
    #center=array['phi'][ny[i]/2,:]

    
    arr=genfromtxt("Sehgal/ca00113and007.csv",delimiter=',')
    
    widths=[]
    velocities=[]
    fig=pylab.figure()
    pylab.plot(arr[:,0]/16.0,arr[:,1]/2.0,'b-.',markersize=9,linewidth=3)
    pylab.plot(arr[:,0]/16.0,arr[:,2]/2.0,'r--',markersize=9,linewidth=3)
    #print len(arr[:,0])
    #print len(smoothList(arr[:,1]))    
    #pylab.plot(arr[:-4,0],smoothList(arr[:,1]),linewidth=3)    
    pylab.ylim(ymin=0.0,ymax=0.5)
    #pylab.xlim(xmax=5.1)

    leg=pylab.legend(labels)
    legtext = leg.get_texts() # all the text.Text instance in the legend
    for text in legtext:
        text.set_fontsize(20) 
    #set(ltext, fontsize='large') # the legend text fontsize

    fig.subplots_adjust(left=0.15,bottom=0.15)  
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
    pylab.xlabel(r'''$x$''',fontsize=30)
    pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(["Simulations","Giavedoni","Heil"],loc=4)
    #pylab.xlim(xmax=15)
    pylab.savefig("bubble_sehgal_new.eps",format="EPS",dpi=300)
    
    
    
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

def Analyze_Vectors():
    from PyNGL import Ngl
    from numpy import genfromtxt
    print os.getcwd()
    capillary_theor=numpy.array([0.03,0.05,0.08,0.1,0.2,0.4,0.6,0.8])
    capillary_str=numpy.array(["3","5","8","10","20","40","60","80"])
 
    good=[0,1,2,3,4,5,6,7]
    good_value=[7]
    style=["b-.","r--"]
    #labels=[r'''$Ca='''+ca+r'''$''' for ca in numpy.array(capillary_theor[good_value],dtype=str)]
    labels=[]    
    #color=["b","r","c","y","b","g"]
    #style_diff=["b-","r:","c-.","y--","b^","g<"]
    
    #ux=array['v'][0]
    #center=array['phi'][ny[i]/2,:]

    for counter,value in enumerate(good_value):
        dir_temp="Results/"+capillary_str[value]
        os.chdir(dir_temp)
        name="capillary200000.npz"
        array=numpy.load(name)
        thicknesses=[]
        dims=array['phi'].shape        
        ux=array['v'][0]
        uy=array['v'][1]
        phase=array['phi']
        x,y=numpy.mgrid[0:dims[0],0:dims[1]]        

        
        center=array['phi'][dims[0]/2,:]
           
        z1 = numpy.min(numpy.where(center < 0.0))
        z2 = numpy.max(numpy.where(center < 0.0))
        if z1==0:
            z2=numpy.min(numpy.where(center>0.0))+dims[0]
            z1=numpy.max(numpy.where(center>0.0))
        print z1,z2
        print "Capillary=",ux[dims[0]/2,z2%dims[1]]*2.0/3.0/math.sqrt(8.0*0.04*0.04/9.0)
        positive=numpy.where(phase>0.0)
        negative=numpy.where(phase<0.0)
        large=numpy.where(numpy.logical_or(x<20,x>dims[0]-20))
        #bounds=numpy.where(numpy.logical_and(y>z1-5,y<z2%dims[2]+10))
        print z1,z2 #,bounds
        #print y<z1-30
        #x_short=x[::3,::25]*deltay
        #y_short=y[::3,::25]*deltax
        #[numpy.where(phase_numpy>0.0)]
        ux=ux-ux[dims[0]/2,z2%dims[1]]        
        #ux[negative]=None
        #ux[large]=None
        #vz_diff_mask[bounds]=None
        
        ux_mask=ux[::5,::100]
        uy_mask=uy[::5,::100]
        x_mask=x[::5,::100]
        y_mask=y[::5,::100]

        
        #pylab.figure(figsize=(20,10))
        #pylab.quiver(y_mask,x_mask,ux_mask,uy_mask,headwidth=6,minlength=0.1)
        #pylab.contour(array['phi'],[0.0],linewidths=[3])
        
        wks_type = "eps"
        wks = Ngl.open_wks(wks_type,"test")
        resources = Ngl.Resources()
     
        #uvar = file.variables["U_GRD_6_ISBL"]
        #vvar = file.variables["V_GRD_6_ISBL"]
        #if hasattr(uvar,"units"):
        #  resources.tiMainString = "GRD_6_ISBL (u,v " + uvar.units + ")"
        #else:
        #resources.tiMainString = "GRD_6_ISBL"
        #if hasattr(uvar,"_FillValue"):
        #    resources.vfMissingUValueV = uvar._FillValue
        # if hasattr(vvar,"_FillValue"):
        #    resources.vfMissingVValueV = vvar._FillValue
        
        
        #resources.tiMainFont    = "Times-Roman"
        #resources.tiMainOn=True
        #resources.tiMainString="Ca=0.22 "
           
        #resources.tiXAxisString = "streamlines"
        resources.vpHeightF = 0.25 # Define height, width, and location of plot.
        resources.vpWidthF  = 3*0.25
        resources.wkPaperSize="A5"
        resources.nglFrame = False
        resources.vfXArray=numpy.linspace(0.0,15.0,len(ux[1,::50]))
        resources.vfYArray=numpy.linspace(0.0,1.0,len(ux[::5,1]))
        #resources.


        
        
        resources2=Ngl.Resources()
        #resources2.tiMainFont    = "Times-Roman"
        #resources2.tiXAxisString = "streamlines"
        #resources2.tiMainOn=True
        #resources2.tiMainString="Ca=0.22"
        
        resources2.wkPaperSize="A5"
        resources2.vpHeightF = 0.25 # Define height, width, and location of plot.
        resources2.vpWidthF  = 3*0.25
        resources2.nglFrame = False
        
        resources2.cnLineLabelsOn = False   # Turn off contour line labels.
        #resources2.cnLinesOn      = False   # Turn off contour lines.
        resources2.cnFillOn       = False    # Turn on contour fill.
        resources2.cnInfoLabelOn   = False 
  
        
        resources2.cnLevelSelectionMode = "ExplicitLevels"  # Select contour levels. 
        resources2.cnMinLevelValF       = 0.0
        #resources2.cnMaxLevelValF       = 0.001
        #resources2.cnLevelSpacingF      = 0.0
        resources2.cnLevelCount=1
        resources2.cnLevels=[0.0]
        #resources2.cnLineThicknesses=[3]
        resources2.cnMonoLineThickness=True
        resources2.cnLineThicknessF=3.0
        
        resources2.lbLabelBarOn=False
        resources2.lbLabelsOn=False
        resources2.sfXArray=numpy.linspace(0.0,15.0,len(ux[1,:]))
        resources2.sfYArray=numpy.linspace(0.0,1.0,len(ux[:,1]))
        
        #plot = Ngl.streamline(wks,uvar[0,::2,::2],vvar[0,::2,::2],resources) 
        #print vz_diff.shape
        #print vy.shape
        x,y=numpy.mgrid[0:dims[0],0:dims[1]]
        #vx=numpy.sin(x)*numpy.sin(y)
        #vy=numpy.cos(x)*numpy.cos(y)
        plot=Ngl.streamline(wks,ux[::5,::50],uy[::5,::50],resources)
        #Ngl.contour(wks,phase[::5,::50],resources2)        
        Ngl.contour(wks,phase,resources2)        
        #plot=Ngl.streamline(wks,vx,vy,resources)
        Ngl.end()
        
        print "Ya zdes'!"
        #print thicknesses[((z1+z2)/2)%dims[1]-z1]        
        #print numpy.std(thicknesses[300:-300])/thicknesses[((z1+z2)/2)%dims[1]-z1]        
        #delta_x=15.0/(len(thicknesses)-1)
        #x=delta_x*numpy.arange(0,len(thicknesses))*float(len(thicknesses))/3000.0
        #pylab.plot(x,thicknesses,style[counter],linewidth=3)
        #labels.append(r'''$Ca='''+str(ux[dims[0]/2,z2%dims[1]]*2.0/3.0/math.sqrt(8.0*0.04*0.04/9.0))[0:4]+r'''$''')
        os.chdir("../..")
    
    #pylab.ylim(ymin=0.0,ymax=0.5)
    #pylab.xlim(xmax=5.1)

    #leg=pylab.legend(labels)
    #legtext = leg.get_texts() # all the text.Text instance in the legend
    #for text in legtext:
    #    text.set_fontsize(20) 
    ##set(ltext, fontsize='large') # the legend text fontsize

    #fig.subplots_adjust(left=0.15,bottom=0.15)  
    #pylab.xticks(fontsize=20)
    #pylab.yticks(fontsize=20)
    #pylab.xlabel(r'''$x$''',fontsize=30)
    #pylab.ylabel(r'''$\delta$''',fontsize=30)
    
    #labels=[r'''$H_{eff}='''+str(value-2)+r'''$''' for value in ny]
    #pylab.legend(["Simulations","Giavedoni","Heil"],loc=4)
    #pylab.xlim(xmax=15)
    #pylab.savefig("bubble_length.eps",format="EPS",dpi=300)

        
        
if __name__=="__main__":
    
    Analyze_Simulations()    
    #Analyze_Velocities()
    #Run_Simulations()
    #Analyze_Bubble()
    #Analyze_Sehgal_Bubble()    
    #Analyze_Vectors()    
    pylab.show()
