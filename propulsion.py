"""
UCSB Rocket Propulsion Laboratory
Propulsion Team Rocket Engine Simulation

This file is converted from MATLAB to python

"""


def calculate_engine_thrust(m_dot, burn_time):
    ########################
    ### Assume pe=1 (accurate if burn is short and majority of combustion occurs at sea level 
    ### Using this data: http://www.braeunig.us/space/comb-OM.htm
    ### Assume chamber pressure, p_c=75 atm (guess for now)
    pe= 101325 #(Pa) (1 atm) exit pressure
    T0=3400 #(K) adiabatic flame temperature @p_c=75atm, mixture ratio = 2.80
    gamma= 1.209 # specific heat ratio @p_c=75atm, mixture ratio = 2.80
    M_bar= 19.7 #(kg/mol) @p_c=75atm, mixture ratio = 2.80

    R_bar= 8314 #(J/(kg*mol*K) universal gas constant
    R= R_bar/M_bar #(J/(kg*K)) specifc gas constant = universal gas constant/ avg molecular weight.


    ### Estimated from pg.4/5 of https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20100040703.pdf
    p0= 1.2*10**6 #(Pa) total pressure of LOX/LCH4, 
    #m_dot= .14 #(kg/s) m_dot=m_ox+m_ch4+m_igniter; #exhaust mass flow rate
    ########################

    Ve= ((2*gamma*R*T0/(gamma-1))*(1-(pe/p0)**((gamma-1)/gamma)))**(1/2) #exit velocity

    T=m_dot*Ve #thrust assuming pe=pinf
    #print('Thrust = ' + str(T) + ' N')
    Isp=T/(9.81*m_dot)
    #print('Isp = ' + str(Isp) + ' s')
    return T
