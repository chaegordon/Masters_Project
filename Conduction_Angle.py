# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:59:33 2019

@author: chaeg
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:59:47 2019

@author: chaeg
"""

"""Slab contour including Conduction, will be slower?"""



"""program to calculate the temp. of paper slab with time varying flux on one surface
and a heat sink on the other"""

#need to check the expression for flux

from scipy.special import erf
import scipy.integrate as integrate
import scipy.special as special
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


"""Can ignore the back flux from paper as will also have forward flux from the
other paper cancelling this"""

l = 4 #micrometers
k = 0.05

#kinematic parameters
g=9.81 #SI
h=0.2
M = 0.5

V = (2*g*h)**(0.5)

#thermal diffusivity
kappa = 8.2*10**-8


#588 #if at 2.8 burst strength
Temp_max = 700 #7300
t_c = 5*10**-6

E_m = 200*10**9
R = 0.05
e = 0.2**0.5


N = np.linspace(0,4,5,dtype = float)
x=np.linspace(0*l,l,10**2)
time=np.logspace(-1,0,10**3, dtype = float)
X,T = np.meshgrid(x,time)


"""ierfc is a function defined in appendix 2 of carslaw and jaeger"""

def erfc(x):
    return 1- erf(x)

def ierfc(x):
    WE= np.exp(-x**2/np.sqrt(np.pi))
    COOL = x*erfc(x)
    return WE - COOL

def dE_1(t):
    counter = np.zeros(T.shape)
    for n in N:
      A = ((2*n+1)*l - X)*(10**-6)/(2*(kappa**(0.5)))
      a1 = ((A/(t*10**-6))**(2))
      a7 = np.exp(-((A)**2)/(t*10**-6))
      a2 = (A/(2*(t*10**-6)**(1.5)))
      a8 = erfc(A*(t*10**-6)**-0.5)
      a5 = ((A/(t*10**-6))**(2))
      a4 = (np.exp(-1*(((A)**(2))/(t*10**-6))))
      a6 = 2*(np.pi)**(-0.5)
      counter += a1*a7 + a2*a8 - a6*a4*a5
    return counter

def dE_2(t):
    counter = np.zeros(T.shape)
    for n in N:
      A = ((2*n+1)*l + X)*(10**-6)/(2*(kappa**(0.5)))
      a1 = ((A/(t*10**-6))**(2))
      a7 = np.exp(-((A)**2)/(t*10**-6))
      a2 = (A/(2*(t*10**-6)**(1.5)))
      a8 = erfc(A*((t*10**-6)**(-0.5)))
      a5 = ((A/(t*10**-6))**(2))
      a4 = (np.exp(-1*(((A)**(2))/(t*10**-6))))
      a6 = 2*(np.pi)**(-0.5)
      counter += ( a1*a7 + a2*a8 - a6*a4*a5)
    return counter

#temperature solution in time independent case per unit flux
def v(t):
    counter = np.zeros(T.shape)
    for n in N:
      A= (2/k)*(kappa*(t*10**-6)**(0.5))
      B=((-1)**(n))
      C=ierfc(((2*n+1)*l-X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      D=-1*ierfc(((2*n+1)*l+X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      counter += ( A*B*(C+D))
    return counter


def d_Theta(t):
    counter =  np.zeros(T.shape)
    for n in N:
      A= (2/k)*(kappa/(t*10**-6))**(0.5)
      B=((-1)**(n))
      C=ierfc(((2*n+1)*l-X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      D=-1*ierfc(((2*n+1)*l+X)*(10**-6)/(2*(kappa*(t*10**-6))**(0.5)))
      A_1 = A*B*(C+D)
      E = (2/k)*(kappa*(t*10**-6))**(0.5)
      F = dE_1(t)
      G = dE_2(t)
      A_2 = E*B*(F + G)
      counter += (A_1 + A_2)
      
    return counter

def collision_time(theta):
    t_c = (1+e)*2.87*(((M/2)**2)/float(R*E_m**2*V*np.cos(theta))**0.2)
    return t_c

def f(t,theta):
    t_c = collision_time(theta)
    flux_time = k*((Temp_max*(np.sin((t*10**-6)*np.pi/t_c))**(0.57)) - d_Theta(T))
    return flux_time


angle = np.linspace(0,1.4,10)

soln = []
max_T =[]

for theta in angle:
        fluxy = np.where(f(T,theta) < k*(588 - d_Theta(T)), f(T,theta), k*(588 - d_Theta(T)))
        #accounts for the time step between points on the convolution, must be 1/len(T)
        delta = 0.001  
        temp = signal.fftconvolve(fluxy,delta*d_Theta(T), mode = 'same')
        #print(f(T).shape)
        #print(d_Theta(T).shape)
        """Theta defined as 0=295K"""
        temperature_re_scaled = temp + 295
        #print(temp.shape)
        soln.append(temperature_re_scaled)
        max_T.append(np.amax(temperature_re_scaled))
        
arrhenius = np.exp(-506/np.array(max_T))
    
angle_deg = np.degrees(angle)


plt.figure()
plt.plot(angle_deg,arrhenius)
plt.title('Arrhenius Factor vs Angle of Incidence')
plt.xlabel('Angle of Incidence / degrees')
plt.ylabel('exp(- 506K / T)')
plt.savefig('Pore_dissipated_Arrhenius_Burst_angle.png')

plt.figure()
CP = plt.contourf(X,T, soln[0], 100)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Void Collapse Heating by Conduction')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')
plt.savefig('Pore_Conduction_Heaviside_Burst.png')

