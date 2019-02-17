# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 18:57:48 2018

@author: Adam

EXIT AREA OFF FROM BOOK BY 4%
"""
import numpy as np
k = 1.1409 #dimensionless cnst ratio of specific heats of gas mixture
pOut = .079154  * 10**6 #0.101325 * 10**6 #pressure at exit of nozzle in mpa
#pChamber = 500 * 6894.76 #np.array([100,200,300,400,500]) * (0.00689476)  * 10**6  #first number is in psi
desiredThrust = 3114 #Newtons
thrustCorrectionFactor = .75
g=9.8
ISP = 300 #seconds
def dims(pChamber,disp):
    pChamber *= 6894.76
    tCoeff= np.sqrt((2*(k**2))/(k-1)*((2/(k+1))**((k+1)/(k-1)))*(1-(pOut/pChamber)**((k-1)/k))) 
    throatArea = desiredThrust/(pChamber*thrustCorrectionFactor*tCoeff) #cm^2 is unit here
    exitAreaRatio = 1/((((k+1)/2)**(1/(k-1)))*((pOut/pChamber)**(1/k))*np.sqrt(((k+1)/(k-1))*(1-((pOut/pChamber)**((k-1)/k)))))
    minChamberArea = 4 * throatArea #book says 3 for appreciable pressure drop, 4 for appreciable chamber velocity
    mDot = desiredThrust/(ISP*g)
    if(False):
        print("expansion ratio is ", exitAreaRatio)
        print("throat area is ", throatArea*10000, "cm^2 at tcoeff = ", tCoeff)
        print("exit area is", (throatArea*exitAreaRatio)*10000,"cm^2 and chamber area should be at least ",minChamberArea*10000,"cm^2")
        print("Mass flowrate is ", mDot, "kg/s")
        print(2*np.sqrt(throatArea/np.pi)*100, "cm is the throat diameter", 2*np.sqrt(minChamberArea/np.pi)*100, "cm is the min chamber diameter", 2*np.sqrt(throatArea*exitAreaRatio/np.pi)*100,"cm is exit diameter")
    return exitAreaRatio,2*np.sqrt(throatArea/np.pi)*100,2*np.sqrt(throatArea*exitAreaRatio/np.pi)*100
    
#print("AreaRatio", exitAreaRatio)
#print(np.array([600,700,800,900,1000])/14.6959)
