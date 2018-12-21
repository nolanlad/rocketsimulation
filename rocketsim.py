from math import pi,pow
import matplotlib.pyplot as plt

METERS2FEET = 3.28084 #meters to feet conversion
P0 = 101.325*1000 #atmospheric pressure in kpa at sea level
T0 = 300.0 #initial temperature in kelvins
g = -9.80665 #gravitational constant
L = 0.0065 #temperature coefficient
R = 8.31447
M = 0.0289644

height_list = [] #create empty lists to reference later
time_list = []

def main():

    fp = open("F23FJ.txt","r") #pull outside file
    line = fp.readline().split(" ") #create two values for each line
    Ft = float(line[0]) #define values in each line
    F  = float(line[1])

    J = 29.8
    m_rocket = 0.74607 #total mass of the rocket
    m_rate = 0.03/2.2 #mass of fuel over time
    burn_time = 2.2 #total burn time
    cd = 0.25539395736100207 #various constants
    r = 0.03163
    rho = 1.225
    thrust = 45.0

    h,t,v = 0.0, 0.0, 0.0 #starting values

    A = pi*r*r
    m = m_rocket
    dt = 1.0e-4 #0.0001 time increment

    while t < burn_time: #calculations while fuel is burning
        T = T0 - L*h #temperature decreases linearly at altitude increases
        P = P0 * pow( 1 - (L*h)/T0, (-1*g*M)/(R*L))  #find changes in pressure as altitude increases
        rho = (P*M)/(R*T) #define rho constant
        if t > Ft:
            line = fp.readline().split(" ")
            Ft = float(line[0])
            F  = float(line[1])
            print((F + g*m - 0.5*rho*v*v*cd*A)/m) #find resultant acceleration

        thrust = F
        a = (F + g*m - 0.5*rho*v*v*cd*A)/m

        v +=(a*dt) #increment velocity by acceleration
        h += (v*dt) #increment height by velocity
        m-=(dt*m_rate) #decrease mass by burn rate
        t+=dt #increase time by dt

        time_list.append(t) #append each value to time list
        height_list.append(h) #append ach value to height list

    while v >= 0.0: #calculations after fuel has burned and rocket still traveling up
        T = T0 - L*h
        P = P0 * pow( 1 - (L*h)/T0, (-1*g*M)/(R*L))
        rho = (P*M)/(R*T)
        a = (g*m - 0.5*rho*v*v*cd*A)/m #no longer an applied force

        v +=(a*dt) #increment values
        h += (v*dt)
        m-=(dt*m_rate)
        t+=dt

        height_list.append(h)
        time_list.append(t)

    print(h*METERS2FEET)
    plt.plot(time_list, height_list, 'r--')
    plt.show()

main()
