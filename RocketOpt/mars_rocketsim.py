"""

UCSB Rocket Propulsion Laboratory
Flight Dynamics Rocket Simulation
Andrew Zakoor
Nolan McCarthy
Adam Poklemba

"""

from math import pi,pow,exp,sqrt,cos,atan
import matplotlib.pyplot as plt
import sys

args = sys.argv

def main(thrust,innerRadius,fuel_mass,dry_rocket_mass):

    def cd_altitude_calc():

        slope_a = ( 224 - 293 ) / 11000 #kelvin decrease per meter

        if height < 11000: #calculates the temperature at a given height
            altitude_temperature = mojave_temperature - ( height * meters_to_feet * temp_slope )
            altitude_density = ( 1.225 ) * ( altitude_temperature / 293 )**(-((9.8)/(slope_a*287)+1))
        else: #above 11k meters the temperature is relatively constant
            altitude_temperature = mojave_temperature - ( 11000 * meters_to_feet * temp_slope )
            pressure_1 = ( 1.013e5 ) * (( altitude_temperature / 293 )**(-9.8/(slope_a*287)))
            density_1 = ( 1.225 ) * ( altitude_temperature / 293 )**(-((9.8)/(slope_a*287)+1))
            pressure_2 = ( pressure_1 ) * real_e**(-(9.8*(height - 11000))/(287*224))
            altitude_density = density_1 * ( pressure_2 / pressure_1 )

        #altitude_pressure = 101.29 * ( altitude_temperature / 288.08 )**5.256 #pressure from temperature
        #altitude_density = altitude_pressure / ( .2869 * altitude_temperature ) #density from pressure
        altitude_mach = 331 + ( 0.6 * ( altitude_temperature - 273.15 ) ) #mach from temperature

        # current_cd = 0.000386 * velocity + 0.334181 #data collected experimentally from SolidWorks, regression line fit to data
        # current_cd = 0.51
        if velocity <= 343:
            current_cd = 0.5 + velocity*0.0004373
        else:
            current_cd = 0.65 * 2.71**(-(-343+velocity)/600)
        return altitude_density, current_cd

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
        force_gravity_list.append(force_gravity)
        force_drag_list.append(force_drag)
        density_list.append(altitude_density)

    def plot_plots():

        plt.subplot(3,2,1)
        plt.plot(time_list, height_list, 'r')
        plt.ylabel('Height (m)')

        plt.subplot(3,2,3)
        plt.plot(time_list, velocity_list, 'b')
        plt.ylabel('Velocity (m/s)')

        plt.subplot(3,2,5)
        plt.plot(time_list, acceleration_list, 'g')
        plt.ylabel('Acceleration (m/s^2)')
        plt.xlabel('Time (s)')

        plt.subplot(3,2,2)
        plt.plot(velocity_list, drag_coefficient_list, 'y')
        plt.ylabel('Drag Coefficient')
        plt.xlabel('Velocity (m/s)')

        plt.subplot(3,2,4)
        plt.plot(time_list, force_gravity_list, 'r--')
        plt.ylabel('Force of Gravity (N)')
        plt.xlabel('Time (s)')

        plt.subplot(3,2,6)
        plt.plot(time_list, force_drag_list, 'b--')
        plt.ylabel('Force of Drag (N)')
        plt.xlabel('Time (s)')

        plt.show()

    """ Fundamental Constants """

    gravitational_constant = 6.67408e-11 #constant capital G
    boltzmann_constant = 1.38064852e-23 #in units of m^2 kg s^-2 K^-1
    gas_constant = 287.05 #J/(kg*K)
    meters_to_feet = 3.28084 #meter to feet conversion
    lbf_to_n = 4.4482216152605 #pound force to newton conversion
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

    """ Rocket Constants """
    #thrust in newtons
    burn_time = 9208 / ( thrust / lbf_to_n ) #estimated total burn time of liquid fuel
    total_rocket_mass = fuel_mass + dry_rocket_mass #total mass
    mass_change = ( fuel_mass / burn_time ) * dt #assuming constant change of mass

    rocket_radius = ( innerRadius / 12 ) * ( 1 / meters_to_feet ) #radius of the rocket, meters
    rocket_area = pi * ( rocket_radius**2 ) #area from radius

    """ Initialize Variables """

    time = 0.0 #sets all variables to zero for start of the flight
    height = 0.0
    velocity = 0.0

    height_track = 0.0
    velocity_track = 0.0

    acceleration_list = [] #sets empty lists to be filled with data
    velocity_list = []
    height_list = []
    time_list = []

    drag_coefficient_list = []
    force_gravity_list = []
    force_drag_list = []
    density_list = []

    while time <= burn_time:

        if time >= 1.0:
            weather_adjust_factor = cos(atan(wind_speed/velocity))

        altitude_density, current_cd = cd_altitude_calc() #finds density and drag coefficient at any given iteration
        force_gravity, force_drag = force_calc() #finds gravity and drag at any time

        if time < 1.0:
            acceleration = ( thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
        else:
            acceleration = ( weather_adjust_factor * thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
            if thrust < ( force_drag + force_gravity ):
                force_drag = thrust - force_gravity
                acceleration = ( weather_adjust_factor * thrust - force_gravity - force_drag ) / total_rocket_mass #finds resultant acceleration
                print("this just ran")

        velocity += ( acceleration * dt ) #increment velocity by small steps
        height += ( velocity * dt ) #same process for height
        time += dt #increase total time by dt
        total_rocket_mass -= mass_change #at this time the rocket is losing fuel mass

        if velocity > velocity_track: #keeps track of the greatest total speed
            velocity_track = velocity

        append_lists() #keep track of all data

    while time > burn_time:

        weather_adjust_factor = cos(atan(wind_speed/velocity))

        altitude_density, current_cd = cd_altitude_calc()
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

    # print("Apogee = ", round(height_track,3), "meters, or ", round(height_track*meters_to_feet,3), "feet.")
    #print("Max velocity was ", round(velocity_track,3), "meters per second.")
    #plot_plots()
    return round(height_track*meters_to_feet,3)

# apogee = main(4400,3.5)
# print(apogee)
