"""

UCSB Rocket Propulsion Laboratory
Flight Dynamics Rocket Simulation
Andrew Zakoor
Nolan McCarthy
Adam Poklemba

"""

from math import pi,pow,exp,sqrt,cos,atan
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import numpy as np
import random
import sys

args = sys.argv

def main(thrust_lbf,innerRadius_in,fuel_mass_lb,dry_rocket_mass_lb):

    def cd_altitude_calc():

        slope_a = ( 224 - 293 ) / 11000 #kelvin decrease per meter

        if height < 11000: #calculates the temperature at a given height
            altitude_temperature = mojave_temperature - ( height * meters_to_feet * temp_slope )
            altitude_density = ( 1.225 ) * ( altitude_temperature / 293 )**(-((9.8)/(slope_a*287)+1))
        else: #above 11k meters the temperature is relatively constant
            altitude_temperature = mojave_temperature - ( 11000 * meters_to_feet * temp_slope )
            pressure_1 = ( 1.013e5 ) * (( altitude_temperature / 293 )**(-9.8/(slope_a*287)))
            density_1 = ( 1.225 ) * ( altitude_temperature / 293 )**(-((9.8)/(slope_a*287)+1))
            pressure_2 = ( pressure_1 ) * real_e**(-9.8*(height - 11000)/(287*224))
            altitude_density = density_1 * ( pressure_2 / pressure_1 )

        altitude_mach = sqrt( gamma * gas_constant * altitude_temperature ) #mach from temperature, in m/s

        cd_v3_list = np.array([0.55432, 0.53437, 0.51071, 0.48687, 0.49139, 0.48243, 0.47812, 0.47577, 0.46904, 0.49098, 0.51463, 0.68300, 0.67543, 0.62601, 0.52901, 0.46041, 0.37132, 0.35000])
        drag_v3_list = np.array([0, 8.77484, 33.31318, 72.91437, 127.70369, 197.82749, 283.87060, 387.60258, 511.45387, 579.75506, 691.77405, 1084.06154, 1267.29417, 1426.57258, 1969.45133, 2425.65958, 3757.03066])
        drag_v4_list = np.array([0, 8.52192, 32.28034, 70.80437, 123.82969, 191.56966, 274.50500, 374.45160, 493.56528, 671.10050, 833.52275, 1338.57703, 1416.65956, 1556.88359, 2128.89588, 2777.19525, 3783.07444])
        v_list = np.array([0, 33, 66, 99, 132, 165, 198, 231, 264, 280, 310, 345, 360, 400, 500, 600, 800, 1200])

        cd_list = []
        for x in range(0,len(v_list)-1,1):
            cd = ( drag_v4_list[x] ) / ( 0.5 * rho_0 * pi * 3.25**2 * 0.0254**2 * ( v_list[x] + 0.001 )**2 )
            cd_list.append(cd)

        mach_list = v_list / solidworks_sos
        re_list = ( v_list * rocket_length * rho_0 ) / solidworks_mu
        current_mach = velocity / altitude_mach #mach at any given point, ranges from 0 - 2.3ish

        if current_mach < 0.3: #defined as the mach at which the mach number plays more role than re
            altitude_mu = mu_0*(altitude_temperature/t_0)**1.5*(t_0+s_0)/(altitude_temperature+s_0)
            current_re = velocity * altitude_density * rocket_length / altitude_mu
            for i in range(0,len(re_list),1):
                if current_re >= re_list[i] and current_re <= re_list[i+1]:
                    m = ( cd_list[i+1] - cd_list[i] ) / ( re_list[i+1] - re_list[i] )
                    current_cd = ( current_re - re_list[i] ) * m + cd_list[i]
                    break
                else:
                    continue
        else:
            for i in range(0,len(mach_list),1):
                if current_mach >= mach_list[i] and current_mach <= mach_list[i+1]:
                    m = ( cd_list[i+1] - cd_list[i] ) / ( mach_list[i+1] - mach_list[i] )
                    current_cd = ( current_mach - mach_list[i] ) * m + cd_list[i]
                    break
                else:
                    continue

        return altitude_temperature, altitude_density, current_mach, current_cd

    def force_calc():

        force_gravity = ( ( gravitational_constant * (total_rocket_mass) * (earth_mass) ) / (mojave_radius + height)**2 )
        force_drag = ( 0.5 * current_cd * altitude_density * velocity**2 * rocket_area )

        return force_gravity, force_drag

    def append_lists():

        acceleration_list.append(acceleration)
        velocity_list.append(velocity)
        height_list.append(height)
        time_list.append(time)

        drag_coefficient_list.append(current_cd)
        current_mach_list.append(current_mach)
        force_gravity_list.append(force_gravity)
        force_drag_list.append(force_drag)
        density_list.append(altitude_density)

    def plot_plots():

        color_list = []

        for x in range(0,6,1):
            color = "%04x" % random.randint(0, 0xFFFF)
            color_2 = '#F' + color + 'F'
            color_list.append(color_2)

        ft_list = np.asarray(height_list) * meters_to_feet
        plt.subplot(3,2,1)
        plt.plot(time_list, ft_list, color_list[0])
        plt.ylabel('Height (ft)')
        plt.suptitle('Insert Rocket Name',fontsize=16)

        fts_list = np.asarray(velocity_list) * meters_to_feet
        plt.subplot(3,2,3)
        plt.plot(time_list, fts_list, color_list[1])
        plt.ylabel('Velocity (ft/s)')

        ftss_list = np.asarray(acceleration_list) * meters_to_feet
        plt.subplot(3,2,5)
        plt.plot(time_list, ftss_list, color_list[2])
        plt.ylabel('Acceleration (ft/s^2)')
        plt.xlabel('Time (s)')

        plt.subplot(3,2,2)
        plt.plot(current_mach_list, drag_coefficient_list, color_list[3])
        plt.ylabel('Drag Coefficient')
        plt.xlabel('Altitude Mach')

        g_lbf_list = np.asarray(force_gravity_list) / lbf_to_n
        plt.subplot(3,2,4)
        plt.plot(time_list, g_lbf_list, color_list[4])
        plt.ylabel('Force of Gravity (lbf)')
        plt.xlabel('Time (s)')

        d_lbf_list = np.asarray(force_drag_list) / lbf_to_n
        plt.subplot(3,2,6)
        plt.plot(time_list, d_lbf_list, color_list[5])
        plt.ylabel('Force of Drag (lbf)')
        plt.xlabel('Time (s)')

        plt.subplots_adjust(left=0.2,wspace=0.4,hspace=0.5,top=0.9)

        plt.show()

    """ Fundamental Constants """

    gravitational_constant = 6.67408e-11 #constant capital G
    boltzmann_constant = 1.38064852e-23 #in units of m^2 kg s^-2 K^-1
    gamma = 1.4
    gas_constant = 287.05 #J/(kg*K)
    solidworks_sos = 345.45 #calculated with values, m/s
    solidworks_mu = 1.8213e-5 #SolidWorks viscosity at 297K
    mu_0 = 1.716e-5 #kg/m-s, standard values`
    t_0 = 273.11 #K, used for altitude_viscosity calc
    s_0 = 110.56 #K, used for altitude_viscosity calc
    rho_0 = 1.1826 #kg/m^3, from SolidWorks
    meters_to_feet = 3.28084 #meter to feet conversion
    lbf_to_n = 4.4482216152605 #pound force to newton conversion
    lb_to_kg = 0.453592
    pa_to_psi = 0.000145038
    real_e = 2.71 #approximate value for e
    sealevel_pressure = 101325 #atmospheric pressure at sea level in Pascals
    earth_mass = 5.972e24 #kilograms
    temp_slope = 2 / 1000 #temperature decreases by 2 degrees for every increase in 1000ft altitude
    dt = 0.05 #increments of time value

    """ Mojave Desert Specific Values """

    mojave_radius = 6371.13738e3 #distance in meters from mojave latitude (35 degrees) to center of earth
    mojave_temperature = 23.8889 + 273.15 #degrees kelvin, approximate for May 2020
    mojave_pressure = 100846.66 #pressure in Pascals, converted from mmHg
    avg_airmass = 28.95 / 6.022e20 #average mass of a single air molecule in kg
    wind_speed = 15. #wind speed in m/s
    rail_height = 60 / meters_to_feet #rail height in meters

    """ Rocket Constants """
    #converts all inputs from imperial to metric
    thrust = thrust_lbf * lbf_to_n #convets thrust to N for calculations
    fuel_mass = fuel_mass_lb * lb_to_kg
    dry_rocket_mass = dry_rocket_mass_lb * lb_to_kg

    burn_time = 9208 / ( thrust_lbf ) #estimated total burn time of liquid fuel, in seconds
    total_rocket_mass = fuel_mass + dry_rocket_mass #total mass, in kg
    mass_change = ( fuel_mass / burn_time ) * dt #assuming constant change of mass

    rocket_radius = ( innerRadius_in / 12 ) * ( 1 / meters_to_feet ) #radius of the rocket, meters
    rocket_area = pi * ( rocket_radius**2 ) #area from radius
    rocket_length = 12.3 / meters_to_feet #in meters
    rocket_roughness = 3e-6 #surface roughness of carbon fiber in meters

    """ Initialize Variables """

    time = 0.0 #sets all variables to zero for start of the flight
    height = 0.0
    velocity = 0.0

    height_track = 0.0
    velocity_track = 0.0
    mach_track = 0.0
    q_track = 0.0

    acceleration_list = [] #sets empty lists to be filled with data
    velocity_list = []
    height_list = []
    time_list = []

    drag_coefficient_list = []
    current_mach_list = []
    force_gravity_list = []
    force_drag_list = []
    density_list = []

    while time <= burn_time:

        if height > rail_height:
            weather_adjust_factor = cos(atan(wind_speed/velocity))

        altitude_temperature, altitude_density, current_mach, current_cd = cd_altitude_calc() #finds density and drag coefficient at any given iteration
        force_gravity, force_drag = force_calc() #finds gravity and drag at any time

        if height <= rail_height:
            acceleration = ( thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
        else:
            acceleration = ( weather_adjust_factor * thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
            if thrust < ( force_drag + force_gravity ):
                force_drag = thrust - force_gravity
                acceleration = ( weather_adjust_factor * thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
                print("Uh oh. Rocket's going down")

        velocity += ( acceleration * dt ) #increment velocity by small steps
        height += ( velocity * dt ) #same process for height
        time += dt #increase total time by dt
        total_rocket_mass -= mass_change #at this time the rocket is losing fuel mass

        if velocity > velocity_track: #keeps track of the greatest total speed
            altitude_mach = sqrt( gamma * gas_constant * altitude_temperature ) #mach from temperature
            velocity_track = velocity
            mach_track = velocity / altitude_mach

        q = velocity**2 * altitude_density * 0.5
        if q > q_track:
            q_track = q

        append_lists() #keep track of all data

    while time > burn_time:

        weather_adjust_factor = cos(atan(wind_speed/velocity))

        altitude_temperature, altitude_density, current_mach, current_cd = cd_altitude_calc()
        force_gravity, force_drag = force_calc()

        acceleration = ( - force_gravity - force_drag ) / total_rocket_mass
        velocity += ( acceleration * dt )
        height += ( velocity * dt )
        time += dt

        height_adjust = weather_adjust_factor * height

        if height > height_track: #stops the simulation when rocket is at apogee
            height_track = height
        else:
            break

        append_lists()

    print('')
    print("Max velocity =", round(velocity_track*meters_to_feet,3), "ft/s")
    print("Max mach =", round(mach_track,3))
    print("Max dynamic pressure =", round(q_track*pa_to_psi,3), "psi")
    print("Max drag force =", round(max(force_drag_list)/lbf_to_n,3), "lbf")
    print("Max acceleration =", round(max(acceleration_list)*meters_to_feet,3), "ft/s^2")
    print("Apogee =" , round(height_track*meters_to_feet,3), "ft")
    print('')
    plot_plots()

thrust_lbf = float(args[1])
rocket_radius_in = float(args[2])
fuel_mass_lb = float(args[3])
dry_rocket_mass_lb = float(args[4])

main(thrust_lbf,rocket_radius_in,fuel_mass_lb,dry_rocket_mass_lb)
