# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 15:58:53 2019

@author: chaeg
"""
import numpy as np
from scipy import optimize
from scipy import integrate

#material parameters
dynamic_viscosity = 18.37*10**-6 #Pas
density = 1.225 #kgm**-3
g = 9.81 #kgms**-1
r = 0.06 #m
density_m = 7810 # kgm**-3
m = density_m*((4/3)*np.pi*r**3)
h = 1.0 #m

#from eqn of motion
t_fin = (2*h/g)**0.5#s

N = 10**3
guess = [0]
time = np.linspace(0,t_fin,N)


def func(x):
    A = 4.91*(x*dynamic_viscosity/(g*(t+10**-8)*density))**0.5
    B = -1*r
    C = (r**2 - x**2)**0.5
    return A + B - C

x_0 = []

for t in time:
  t = t
  A = optimize.fsolve(func, x0=guess)
  x_0.append(A)

def part_int(t):
    B = (np.pi*((x_0[int((t/t_fin))])**2)*density*(g*t)**2)/(2*m)
    return B

#solving coupled first order ODEs in place of 2nd order ODE

def dy_dt(U, t):
    # Here U is a vector such that y=U[0] and z=U[1]. This function should return [y', z']
    return [U[1], part_int(t)]
U0 = [0, 0]

Us = integrate.odeint(dy_dt, U0, t=time)
ys = Us[:,0]

venturi_y = ys[N-1]
print('The ball drifts ' + str(venturi_y) + 'm')


