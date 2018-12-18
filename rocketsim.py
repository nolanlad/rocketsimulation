from math import pi,pow

METERS2FEET = 3.28084
P0 = 101.325*1000
T0 = 300.0
g = -9.80665
L = 0.0065  
R = 8.31447
M = 0.0289644

def main():

    fp = open("F23FJ.txt","r")
    line = fp.readline().split(" ")
    Ft = float(line[0])
    F  = float(line[1])

    J = 29.8
    m_rocket = 0.74607
    m_rate = 0.03/2.2
    burn_time = 2.2
    cd = 0.25539395736100207
    r = 0.03163
    rho = 1.225
    thrust = 45.0

    h,t,v = 0.0, 0.0, 0.0

    A = pi*r*r
    m = m_rocket
    dt = 1.0e-4
    while t < burn_time:
        T = T0 - L*h
        P = P0 * pow( 1 - (L*h)/T0, (-1*g*M)/(R*L))
        rho = (P*M)/(R*T)
        if t > Ft:
            line = fp.readline().split(" ")
            Ft = float(line[0])
            F  = float(line[1])
            print((F + g*m - 0.5*rho*v*v*cd*A)/m)
            
        thrust = F
        a = (F + g*m - 0.5*rho*v*v*cd*A)/m
        
        v +=(a*dt)
        h += (v*dt)
        m-=(dt*m_rate)
        t+=dt

    while v >= 0.0:
        T = T0 - L*h
        P = P0 * pow( 1 - (L*h)/T0, (-1*g*M)/(R*L))
        rho = (P*M)/(R*T)
        a = (g*m - 0.5*rho*v*v*cd*A)/m
        
        v +=(a*dt)
        h += (v*dt)
        m-=(dt*m_rate)
        t+=dt

    print(h*METERS2FEET)
    

main()



