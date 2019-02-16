# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 20:23:01 2019

@author: Adam
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import PcOpt as alt
from mpl_toolkits.mplot3d import Axes3D


def f(pc,m,d):
    return -alt.pcOpt(pc,m,d)

count = 10
mass = np.linspace(10,75,count)
diameter = np.linspace(5.5,8,count)
allowedMass = []
pc =  np.zeros((count,count))
altitude = np.zeros((count,count))



for di,d in enumerate(diameter):
    maxmass = True
    for mi,m in enumerate(mass):
        cP = opt.fminbound(f,150,600,args=(m,d),xtol = 1)
        pc[di][mi] = cP
        altitude[di][mi] = alt.pcOpt(cP,m,d)     
        print(mi+(di*count))
        if(altitude[di][mi] < 45000 and maxmass):
            allowedMass.append([d,m]) #alt.pcOpt(cP,i,d,moreData=True)[1]])
            maxmass = False

am = np.array(allowedMass)
x,y = np.meshgrid(mass,diameter)

plt.figure(1)
plt.plot(am[:,0],am[:,1])
plt.xlabel("Rocket Diameter(in)")
plt.ylabel("Dry Mass (lb)")

plt.figure(2)
ax = plt.axes(projection='3d')
ax.set_xlabel("Non-Propulsion Mass (lb)")
ax.set_ylabel("Diameter (in)")
ax.set_zlabel("Optimal Pc (psi))")
ax.plot_surface(x,y, pc, cmap='viridis', edgecolor='none')



plt.figure(3)
ax = plt.axes(projection='3d')
ax.set_xlabel("Non-Propulsion Mass (lb)")
ax.set_ylabel("Diameter (in)")
ax.set_zlabel("Apogee (ft))")
ax.plot_surface(x,y, altitude, cmap='viridis', edgecolor='none')

