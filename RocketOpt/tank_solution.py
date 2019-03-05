from math import sqrt, pow, pi
from scipy.optimize import fsolve


class TankSolution(object):

  def __init__(self, tank_pressure, mass_meth, mass_lox,disp=False):

    self.pressure = tank_pressure  # [psi]

    self.r_outer = 6.3/2  # [in] r from centerline out to the inside of the rocket wall

    dens_alum = .098 # lbs/in^3

    dens_lox = 63.91  # [lbs/ft^3] density of lox
    # note that x**y is equivalent to x^y in matlab
    dens_lox = dens_lox / 12 ** 3  # converted to lbs/in^3
    self.V_lox = mass_lox / dens_lox  # [in^3]

    dens_meth = 26.64  # [lbs/ft^3] density of meth
    dens_meth = dens_meth / 12 ** 3  # converted to lbs/in^3
    self.V_meth = mass_meth / dens_meth  # [in^3] the volume of the fuel

    y_stgth_al = 50038.02  # [psi] yield strength of aluminum

    self.t_wall = sqrt(3) / 2 * self.pressure * self.r_outer / y_stgth_al

    self.ro_lox = self.r_outer - self.t_wall  # the radius up to far edge of space inside lox tank

    # solve for the point where tanks meet and the length of the tanks
    initial_guess = (self.r_outer / 2, 1)
    self.xmid, self.L = fsolve(self.system, initial_guess)

    # now calculate the mass of each tank ( surface area * t_wall)
    r_meth = self.xmid - self.t_wall/2
    self.V_tank_shell_meth = 2* (4*pi*r_meth**2*self.t_wall) + 2*pi*r_meth*self.L*self.t_wall
    self.mass_tank_meth = self.V_tank_shell_meth * dens_alum

    ri_lox = self.xmid + self.t_wall/2
    ro_lox = self.r_outer - self.t_wall/2
    V_i_cyl_shell = 2*pi*ri_lox*self.L*self.t_wall
    V_o_cyl_shell = 2*pi*ro_lox*self.L*self.t_wall
    r_tor = self.r_outer - self.t_wall/2 - (self.xmid + self.r_outer)/2
    R_tor = (self.xmid + self.r_outer)/2
    V_tor_shell = 2*pi*r_tor * 2*pi*R_tor * self.t_wall
    self.V_tank_shell_lox = V_i_cyl_shell + V_o_cyl_shell + V_tor_shell
    self.mass_tank_lox = dens_alum * self.V_tank_shell_lox


  # this function defines the system of equations relating xmid, volumes, lengths, etc
  # scipy's fsolve() function will use it to numerically solve the nonlinear system
  # p = [xmid, L_lox]
  def system(self, p):
    xmid, L = p

    ri_meth = xmid - self.t_wall  # radius to the end of empty space inside methane (inner) tank
    ri_lox = xmid + self.t_wall
    r_tor = (self.r_outer - xmid - 2*self.t_wall)/2 # radius of toros circle
    R_tor = (xmid + self.r_outer) / 2  # radius of the sweep of the toroid
    V_tor = (pi * r_tor ** 2) * (2 * pi * R_tor)  # cross sectional area * circumference of sweep

    L_lox = (self.V_lox - V_tor) / (pi * (self.ro_lox ** 2 - (ri_lox) ** 2))

    L_meth = (self.V_meth - 4 / 3 * pi * ri_meth ** 3) / (pi*ri_meth ** 2)

    # both L's must = L
    f1 = L - L_meth
    f2 = L- L_lox


    return (f1, f2)