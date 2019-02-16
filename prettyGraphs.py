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

count = 20
mass = np.linspace(10,75,count)
diameter = np.linspace(5.5,8,count)
allowedMass = []
pc =  np.zeros((count,count))
altitude = np.zeros((count,count))


c = -1
t = -1
for j in diameter:
    c+=1
    t=-1
    maxmass = True
    for i in mass:
        t += 1
        cP = opt.fminbound(f,150,600,args=(i,j),xtol = 1)
        pc[c][t] = cP
        altitude[c][t] = alt.pcOpt(cP,i,j)     
        print(t+(c*count))
        if(altitude[c][t] < 45000 and maxmass):
            allowedMass.append([j,i]) #alt.pcOpt(cP,i,j,moreData=True)[1]])
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

