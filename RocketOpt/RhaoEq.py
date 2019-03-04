# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:33:31 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import NozzleCreator as nozz
import csv

#####Data and Rhao process from http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
def design(pC,disp = False,reso = .01):
    #Data 
    eRat = np.array([3.5,4,5,6,7,8,9,10,20])
    tNData = np.array([19.8,21.6,23,24,24.8,25.3,26,27.2,28.9])
    tEData = np.array([14.8,14,13,12,11.4,11,10.8,10.2,9])
    tNC = np.polyfit(eRat,tNData,2)
    tNF = np.poly1d(tNC)
    tEC = np.polyfit(eRat,tEData,2)
    tEF = np.poly1d(tEC)
    #inputs
    minT = 1
    res = .01
    #nozzle and throat characteristics from my code
    epsilon,dT,dE = nozz.dims(pC,disp)
    rE,rT = dE/2,dT/2
    lN = .8*(rT*(np.sqrt(epsilon)-1))/np.tan(np.deg2rad(15))
    tN,tE = tNF(epsilon),tEF(epsilon)
    
    #all inputs in degrees and cm, outputs in cm
    def entranceT(theta): #converging
        theta = np.deg2rad(theta)
        x = 1.5 * rT * np.cos(theta) 
        y = 1.5 * rT * np.sin(theta) + 2.5 * rT
        return x,y
    
    def exitT(theta): #diverging
        theta = np.deg2rad(theta)
        x = .382 * rT * np.cos(theta)
        y = .382 * rT * np.sin(theta) + 1.382 * rT
        return x,y
    
    def nozzle(t): #nozzle
        nX,nY = exitT(tN-90)
        eX,eY = lN,rE
        m1,m2 = np.tan(np.deg2rad(tN)),np.tan(np.deg2rad(tE))
        c1,c2 = nY - m1*nX, eY - m2*eX
        qX,qY = (c2-c1)/(m1-m2),(m1*c2 - m2*c1)/(m1-m2)
        x = (1-t)**2 * nX + 2*(1-t)*t*qX + t**2*eX
        y = (1-t)**2 * nY + 2*(1-t)*t*qY + t**2*eY  
        return x,y
    #find theta for which converging width = nozzle exit width
    def diff(theta):
        x1,y1 = entranceT(theta)
        x2,y2 = nozzle(1)
        return y1-y2
       
    for i in np.arange(-300,-100,.5):
        if diff(i) < 0 :
            minT =  i
            break
    
    thetaE = np.arange(minT,-90+res,res)
    x,y = entranceT(thetaE)
    
    thetaEx = np.arange(-90, tN - 90+res, res)
    x1,y1 = exitT(thetaEx)
    
    t = np.arange(0,1+res*10,res*10)
    x2,y2 = nozzle(t)
    
    vCon = np.pi * np.trapz(y**2,x)
    vCyl = 100*(rT**2*np.pi) - vCon
    lCh = vCyl/(rE**2*np.pi)
    xC,yC = [x[0]-lCh],[y[0]]
    
    xNet = np.concatenate((xC,x,x1,x2))
    yNet = np.concatenate((yC,y,y1,y2))
    
    xInterp = np.arange(xC[0],x2[-1],reso)
    yInterp = np.interp(xInterp,xNet,yNet)
    def plot():
        plt.figure(1)
        plt.plot(xNet,yNet)
        plt.plot(xInterp,yInterp)
#        plt.plot(xInterp,(yInterp/min(yInterp))**2)
    if(disp):
       plot()
       pass
    return xInterp,yInterp




xVals,yVals = design(550,reso=.1)
xVals *= .01
xVals -= xVals[0]
yVals *= .01

with open('Nozzle550psi.csv',mode='w',newline='') as nozzFile:
    writer = csv.writer(nozzFile)
    writer.writerow(xVals)
    writer.writerow(yVals)

