"""

UCSB Rocket Propulsion Laboratory
Flight Dynamics Rocket Simulation
Andrew Zakoor

"""

from math import pi,pow
import matplotlib.pyplot as plt


""" Fundamental Constants """

gravitational_constant = 6.67408 * ( 10**(-11) ) #constant capital G
meters_to_feet = 3.28084 #meter to feet conversion
lbf_to_n = 4.4482216152605 #pound force to newton conversion
dt = 1 #increments of time value

""" Mojave Desert Specific Values """

mojave_radius = 6371.13738 * ( 10**3 ) #distance in meters from mojave latitude (35 degrees) to center of earth
earth_mass = 5.972 * ( 10**24 ) #kilograms
mojave_temperature = 23.8889 #degrees celsuis, approximate

""" Rocket Constants """

thrust = 1000.0 * lbf_to_n #in Newtons, assuming a approximate and approximately constant force
burn_time = 10.0 #estimated total burn time of liquid fuel

fuel_mass = 15.6 #kilograms of liquid oxygen/methane
dry_rocket_mass = 45 #actual value unknown at the moment
total_rocket_mass = fuel_mass + dry_rocket_mass #total mass
mass_change = ( fuel_mass / burn_time ) * dt #assuming constant change of mass

rocket_radius = 0.1 / 2 #radius of the rocket, meters
drag_coefficient = 0.25 #need to find an actual value
area = pi * ( rocket_radius**2 )

""" Initial Variables """

time = 0.0
height = 0.0
velocity = 0.0

height_list = [] #create empty lists to reference later
time_list = []

while time <= burn_time:

    force_gravity = ( ( gravitational_constant * (total_rocket_mass) * (earth_mass) ) / (mojave_radius + height)**2 )

    # density = pressure / ( gas_constant * temperature )
    # force_drag = ( 0.5 * drag_coefficient * density * velocity**2 * area )

    acceleration = (thrust - force_gravity) / total_rocket_mass

    velocity += ( acceleration * dt )
    height += ( velocity * dt )
    time += dt
    total_rocket_mass -= mass_change

    height_list.append(height)
    time_list.append(time)

while time > burn_time:

    force_gravity = ( ( gravitational_constant * (total_rocket_mass) * (earth_mass) ) / (mojave_radius + height)**2 )

    acceleration = ( - force_gravity) / total_rocket_mass

    velocity += ( acceleration * dt )
    height += ( velocity * dt )
    time += dt

    height_list.append(height)
    time_list.append(time)

    if height <= 0:
        break

plt.plot(time_list, height_list, 'r--')
plt.show()
