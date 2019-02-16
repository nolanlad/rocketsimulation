"""

UCSB Rocket Propulsion Laboratory
Flight Dynamics Rocket Simulation
Nolan McCarthy

"""

""" Geometric Rocket Constants """

from fin import *
from math import pi

rocket_radius = 0.1 / 2 # meters
fin_width     = 0.4  # meters
fin_height    = 0.4  # meters
rocket_height = 4.0     # meters
center_of_gravity = 1.0
center_of_gravity = rocket_height/2

mass = 20.0 #kg
I = 1/12. * (mass*rocket_height**2)

""" Construct Fin Geometry """
fin1  =  Triangle([0,0,0],[0,0,fin_height],[0,fin_width,0])
fin1 = triadd([0,0.1,0],fin1)
rot180 = rot([0,0,1],pi) # rotation matrix for 180 degrees around z-axis
fin2  =  trimatmul(rot180,fin1) # make second fin mirror of first
F = Fin()
F.append(fin1)
F.append(fin2)
F = finadd([0,0,-1*center_of_gravity],F) # move fin so that it rotates around center of gravity


""" initial offset of rocket """
angle = np.pi/24

rot_mat = rot([0,1,0],angle)
F2 = finmatmul(rot_mat,F)
ang = [angle]
alpha = 0.0
omega = 0.0
dt = 0.01

f_accum = 0.0

wind = [2.0,0,0]

for t in np.arange(0,90,dt):
    c, f = fin_drag(F2,[0,0,-50],1.23)
    c_wind, f_wind = fin_drag(F2,[-15,0,0],1.23)
    c = c[0]
    torque = c*f - c_wind[1]*f_wind
    alpha = torque/I
    omega += (alpha*dt)
    angle += (omega*dt)
    rot_mat = rot([0,1,0],angle)
    F2 = finmatmul(rot_mat,F)
    ang.append(angle)


plt.plot(ang)
plt.ylabel("angular displacement (radians)")
plt.show()

