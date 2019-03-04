# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 17:16:22 2019

@author: Adam P, Andrew C
"""
import numpy as np
import mars_rocketsim as v1
import RhaoEq as nProf
import scipy.optimize as opt



#inputs in lb,in
def pcOpt(cPress,mass,diameter,moreData = False):
    nPropMass = mass * 0.453592 #mass of elements of rocket not calculated here
    d = diameter * .0254 #rocket diameter in m
    rHe = 2077.1 #J/kg K
    It = 9208 #Total Impulse (9208 lb-sec Max)
    It = 4.44822*It #Total Impulse (Ns)
    MR = 2.6 #mixture ratio = mdot_o/mdot_f @Pc=25atm #http://www.braeunig.us/space/comb-OM.htm
    
    parameters,values = [],[]
    
    IspD2 = np.array([2169.8,2472.9,2624.9,2721.5,2790.9,2844.2,2887.1,2922.8,2953.2,2979.6]) # ISP CEA data
    pcD2 = np.array([100,200,300,400,500,600,700,800,900,1000])
    
    f2 = np.polyfit(pcD2,IspD2,8)
    q2 = np.poly1d(f2) #polynomial least squares fit for ISP data
    
    #clean input
    pc = 1
    if(type(cPress) != 'list' and type(cPress) != 'numpy.array'):
        pc = np.array([cPress]) #np.array([300]) #
    else:
        pc = cPress
    Isp = q2(pc)

    pd_inj = 0.2*pc #pressure drop across injector 
    pd_cool = 0.2 * pc #pressure drop across cool ing channels
    pd_plumb = 0.2*pc #pressure drop across plumbing
    pt = pc+pd_inj+pd_cool+pd_plumb #pressure in fuel and ox tanks
    pt_si = pt*6894.757 #chamber pressure (Pa)
    pr = pt_si * 1.5 #rated pressure (Pa) --> safety factor of 1.5
    
    rho_al = 2700 #density of Al 6061 T6 (2700 kg/m**3)
    s_al = 290*10**6 #yield stress of Al 6061 T6 at -80C (290 MPa) ...
    s_cu = 70 * 10**6 #at 700k
    rho_cu = 8933 * 10**-6 #kg/m^3 to kg/cm^3
        #http://www.matweb.com/search/datasheet_print.aspx?matguid=1b8c06d0ca7c456694c7777d9e10be5b
    
    def propellantsMassVols(disp = False):
        mp = (It)/Isp #propellant mass
        mf = mp/(1+MR) #fuel mass
        mo = mp/(1+1/MR) #ox mass
        Vo = mo/1141 #ox volume (density of liquid oxygen = 1141 kg/m**3)
        Vf = mf/425.6 #fuel volume (density of liquid methane = 425.6 kg/m**3)
        mp = mo+mf #propellant mass
        if(disp):
            parameters.extend(["Mass of Fuel (kg)","Mass of Oxidiser","Mass of Propellants",
                               "Volume of Fuel (l)","Volume of Oxidiser","Volume of Propellants","Mdot (kg/s)"])
            values.extend([mf,mo,mp,Vf*1000,Vo*1000,(Vf+Vo)*1000,mp/(9208/Isp)])
        return mp,Vo,Vf,mo,mf
    
    def propellantTanksMass(oxVol, methVol, disp = False):
        th_t = (pr*(d/2))/(s_al) #thickness of prop/ox tanks required
        
        A_t = np.pi*((d/2)-th_t)**2 #inner cross sectional area for fuel and ox
        l_ot = oxVol/A_t #length of ox tank
        l_ft = methVol/A_t #length of fuel tank
        
        V_ot = 2*(np.pi*(d/2)**2)*th_t + (np.pi*((d/2)**2-(d/2-th_t)**2))*l_ot #material volume of ox tank
        V_ft = 2*(np.pi*(d/2)**2)*th_t + (np.pi*((d/2)**2-(d/2-th_t)**2))*l_ft #material volume of fuel tank
        
        mot = V_ot*rho_al #ox tank mass
        mft = V_ft*rho_al #fuel tank mass
        mt = mot+mft #total tank mass
        if(disp):
            parameters.extend(["Oxidiser Tank Mass (kg)","Fuel Tank Mass","Total Tank Mass","Tank Thickness(cm)"])
            values.extend([mot,mft,mt,th_t*100])
        return mt
    
    def pressurizingSysMass(oxVol,methVol,disp = False):
        mtHe = []
        mHe = []
        vOut = (4/3)*np.pi*(d/2)**3
        def volumeHe(pressHe):
            return (4/3)*np.pi*((d/2)*(1-pressHe*1.5/(2*s_al)))**3
        def hePressure(press,ind):
            return pr[ind]*((oxVol[ind]+methVol[ind])+volumeHe(press))-press*volumeHe(press)
        for i in range(0,len(pc)):
            pHe = opt.ridder(hePressure,10000,100000000,args=i)
            vHe = volumeHe(pHe)
            mHe.append((vHe*pHe/(rHe*300))*4)
            mtHe.append((vOut-vHe) * rho_al)
        if(disp):
            parameters.extend(["He Volume(l)","He Pressure (MPA)","He Tank Mass(kg)","He Mass"])
            values.extend([vHe,pHe,mtHe,mHe])
        return mtHe + mHe

    def engineMass(disp = False):
        mEng = []
        tEng = []
        for i in pc:
            xD,yD = nProf.design(i,disp = disp) #coordinates from RHAO Code
            vIn = np.pi * np.trapz(yD**2,xD) #inner volume
            engineThickness = (i*6894.757 * 1.5 * max(yD) * 2 )/s_cu
            m = np.array(-1 * np.gradient(yD,xD)) 
            m = np.where(np.isnan(m),0,m)
            xOut,yOut = xD + engineThickness/np.sqrt(1+m**2)*m, yD+engineThickness/np.sqrt(1+m**2)# outside wall coords
            vOut = np.pi * np.trapz((yOut)**2,xOut) #outside volume
            vEng = vOut - vIn 
            mEng.append(vEng * rho_cu)
            tEng.append(engineThickness)
            print("SA", np.pi*np.trapz(yD,xD))
            print(vIn)
        if(disp):
            parameters.extend(["Engine Mass (kg)","Engine Thickness (cm)","Exit Diameter","Throat Diameter","Engine Length","Chamber Length"])
            values.extend([mEng,tEng,2*max(yD),2*min(yD),max(xD)-min(xD),xD[np.argmin(yD)]-xD[0]])
        return mEng
    
    apogee = []
    mp,oxVol,methVol,oxMass,methMass = propellantsMassVols(disp=moreData)
    mt = propellantTanksMass(oxVol,methVol,disp=moreData)
    #mHe = pressurizingSysMass(oxVol,methVol,disp=True)
    mEng = engineMass(disp=moreData)
    for i in range(0,len(pc)):
        apogee.append(v1.main(3114,diameter/2,mp[i],nPropMass+mt[i]+mEng[i]))#+mHe[i]))

    dryMass =  nPropMass+mt[0]+mEng[0] #+mHe[0]
    wetMass = dryMass + mp[0]
    if(not moreData):
        return apogee[0]
    else:
        parameters.extend(["Apogee (ft)","Dry Mass(lb)","Wet Mass","Isp (s)"])
        values.extend([apogee,dryMass*2.2,wetMass*2.2,Isp/9.8])
        print("       Parameter        |     Value     ")
        for parameter, value in list(zip(parameters, values)):
            print(parameter + " "*(24-len(parameter))+"| " + str(value))
        return apogee,dryMass * 2.2,wetMass*2.2,Isp/9.8


pcOpt(300,45,6.5,moreData=True) 


#print(v1.main(3114,3.25,16.3,25))   
#print("Max apogee ", max(apogee), " ft at Pc ",  pc[np.argmax(apogee)], " psi dry mass ", dryMass , " kg wet mass ", wetMass)
#plt.plot(pc,apogee)
#plt.ylabel("Apogee (ft)")
#plt.xlabel("Pc (psi)")

