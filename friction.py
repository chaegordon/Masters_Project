# -*- coding: utf-8 -*-
"""
Created on Thu May  2 16:46:57 2019

@author: chaeg
"""

from scipy.special import erfc
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

"""I think I need to consult a moving source or smth because while work is done 
during slip as the ball moves the place its heating will too! Then can use this as
the same model as for grit and say that its likely single grit is sufficient as 
its melting point high enough and model shows can reach t he temp very quickly"""

#kinematic parameters
g=9.81 #SI
h=0.2
m = 0.5

V = (0.5*g*h*(3)**(0.5))**(0.5)
#material properties
E=200*10**9
k_s = 40.1 #thermal conductivity W/mK
k_p = 0.05
sigma = 2*10**9 #Pa using value for steel
kappa_p = 8.2*10**-8 #m^2/s thermal diffusivity
kappa_m = 1.172*10**-5
b =0.002  # contact radius in m with the PAPER, to be adjusted when calculate more precise 

tau = 400*10**3 # using burst strength of paper
x_0 = 6.75*10**-5 #max comp. steel on steel sphere

Pe_p = (V*b)/(2*kappa_p)
Pe_m = (V*b)/(2*kappa_m)

eta = 1 / (1 + (kappa_m*(Pe_m)**(0.5))/(kappa_p*(Pe_p)**(0.5)))

l = 100 #micrometers is the paper width
N = np.linspace(0,4,5,dtype = float)
x=np.linspace(0*l,l,10**2)
time=np.logspace(-2,-1,10**3, dtype = float) #in microseconds
X,T = np.meshgrid(x,time)

R = 0.05
R_2 = 0.05
M = 0.5
mew = 0.1*(5*10**-5)/(2*b) # this is the x_0/length of paper ie. the deflection of the paper
e = 0.2**0.5 #this is the normal hysteresis COR

#NEED t-C calculated each time
t_c = 5*10**-6
t_f = (1+e)*t_c

#material stiffness constants
E_m = 200*10**9
k_m = (((1-0.295**2)/np.pi*E_m) +((1-0.295**2)/np.pi*E_m))**-1*(4/3*np.pi)*((R*R_2)/(R+R_2))**0.5
E_p = 2*10**9 
G_paper = E_p/2*(1+1/3)
squiggle_2 = k_m/G_paper

angle = np.linspace(0,1.4,10)
w=125600

soln = []
max_T =[]
Dissipated_E = []

for theta in angle:
    t_c = (1+e)*2.87*(((M/2)**2)/float(R*E_m**2*V*np.cos(theta))**0.2)
    Area = np.sin(theta)*2*np.pi*R*V/w + 2*t_c*np.cos(theta)*V**(1.5)*(2*R*np.sin(theta)/w)**0.5
    x_1 = -1*R*np.sin(theta)
    x_2 = R*np.cos(theta)
    #components are for the inertia tensor of THE COLLIDING BALL
    m_1 = M*R**2/(R**2+x_1**2)
    m_2 = M*R**2/(R**2+x_2**2)
     #need to make this E-eff
    t_c = 2.87*(((m_1*m_2/(m_1+m_2))**2)/R*E_m**2*V*np.cos(theta))**0.2
    #natural frequencies of oscillation normal and tangential respectively
    omega=((G_paper/2*m_1)*((1+(m_1*squiggle_2/m_2))+(((1-(m_1*squiggle_2/m_2))**2)+(4*m_1**2/G_paper**2)*x_1**2*x_2**2/M**2*R**4)**0.5))**0.5
    omega_t=((G_paper/2*m_1)*((1+(m_1*squiggle_2/m_2))-(((1-(m_1*squiggle_2/m_2))**2)+(4*m_1**2/G_paper**2)*x_1**2*x_2**2/M**2*R**4)**0.5))**0.5
    #useful constant factor
    a=(omega/omega_t)/(squiggle_2*np.tan(theta)*mew)
    #initial guesses for roots
    x_nought = [0.01, 0.8, 1.5,2.2]
    def find_t1(theta): 
        if np.tan(theta) < mew*squiggle_2:
            def func_t1(i):
                return a*np.sin(omega_t*i) - np.sin((omega/e)*i+(np.pi/2)*(1-1/e))
            t_1 = optimize.fsolve(func_t1,x0=x_nought)
        else:
            def func_t1(i):
                return a*np.sin(omega_t*i) - np.sin(omega*i)
            t_1 = optimize.fsolve(func_t1,x0=x_nought)
        return t_1

    def find_t2(theta):
        if np.tan(theta) < mew*m_2/m_1:
            def func_t2(i):
                return (-1/mew)*np.tan(theta) + (m_2/m_1) - np.cos(omega*i)*(squiggle_2 - m_2/m_1)
            t_2 = optimize.fsolve(func_t2,x0=x_nought)
        else:
            def func_t2(i):
                return (-1/mew)*np.tan(theta) + (m_2/m_1) - np.cos((omega/e)*i+(np.pi/2)*(1-1/e))*(squiggle_2 - m_2/m_1)
            t_2 = optimize.fsolve(func_t2,x0=x_nought)
            return t_2
    g_1 = e*(1 + np.cos((omega*find_t1(theta)[0]/e)+ (np.pi/2)*(1 - 1/e)))
    v_2_0 = np.cos(theta)*V
    v_1_0 = np.sin(theta)*V


    #delete these when can
    v_1_1 = v_1_0*np.cos(omega_t*find_t1(theta)[0])
    v_1_tf = v_1_1 + (m_2/m_1)*g_1*mew*v_2_0
    v_2_tf  = e*v_2_0*np.cos((omega*t_f/e)+ (np.pi/2)*(1 - 1/e))

    #takes argument t which gives which time segment it is referring to

    def v_1_2(theta):
        if find_t2(theta)[0] <= t_c:
            v_1_2 = v_1_0 - mew*(m_1/m_2)*v_2_0*(1-np.cos(omega*find_t2(theta)[0]))
            return v_1_2
        else:
            v_1_2 = v_1_0 - mew*(m_1/m_2)*v_2_0*(1-np.cos(((omega/e)*find_t2(theta)[0])+(np.pi/2)*(1 - 1/e)))
            return v_1_2

    def find_t3(theta): 
        def func_t3(i):
            B_1 = (omega*m_2/(m_1*omega_t**2))
            B_2 = np.sin(omega*find_t2(theta)[0])*np.cos(omega_t*(i-find_t2(theta)[0]))
            B_3 = (omega*v_1_2(theta)/(omega_t*mew*v_2_0))
            B_4 = np.sin(omega_t*(i - find_t2(theta)[0]))
            B_5 = squiggle_2*np.sin((omega/e)*i+(np.pi/2)*(1-1/e))
            A = np.abs(B_1*B_2 - B_3*B_4) - B_5
            return A
        t_3 = optimize.fsolve(func_t3,x0=x_nought)
        return t_3

    def v_1_3(theta):
        A_1 = omega_t* (mew*v_2_0*omega*m_2/(m_1*omega_t**2))
        A_2 = np.sin(omega*find_t2(theta)[0])
        A_3 = np.sin(omega_t*(find_t3(theta)[0]-find_t2(theta)[0]))
        A_4 = v_1_2(theta)*np.cos(omega_t*(find_t3(theta)[0]-find_t2(theta)[0]))
        return A_1*A_2*A_3 + A_4

    def v_1_4(theta):
        if mew*squiggle_2 < np.tan(theta) < mew*(1+e)*(m_2/m_1):
            g_3 = e*(1+np.cos((omega/e)*find_t3(theta)[0]+(np.pi/2)*(1-1/e)))
            v_1_4 = mew*v_2_0*(v_1_3(theta)/(mew*v_2_0) + (m_2/m_1)*g_3)
        else:
            v_1_4 = v_1_0 - mew*v_2_0*(m_2/m_1)*(1+e)
    
        return v_1_4
    

    #Calculates the fraction of energy dissipated by friction

    def D_frac(theta):
        if np.tan(theta) < mew*squiggle_2:
            A= -2*mew**2*g_1*((np.tan(theta)/mew)*np.cos(omega_t*find_t1(theta)[0]) + 0.5*g_1*(m_2/m_1))
            
        elif mew*squiggle_2 < np.tan(theta) < mew*(1+e)*(m_2/m_1):
                
                A= mew**2*(m_1/m_2)*( (np.tan(theta)/mew)**2 - (v_1_2(theta)/mew*v_2_0)**2 + (v_1_3(theta)/mew*v_2_0)**2 -(v_1_4(theta)/mew*v_2_0)**2)
        else:
            A=2*mew**2*(1+e)*((np.tan(theta)/mew - (m_2/2*m_1)*(1+e)))
        return A

    #once this adjusted to be only the tangential lost E will be of corrrect order of mag    
    q_in = e**2 * (m*g*h)*(np.sin(theta))**2*D_frac(theta)/((np.pi*b**2)*10**-5)

    def Temp(t):
    
        A=2*eta*q_in/k_p #heat/t/A into the slab
        B=(kappa_p*t*10**-6)**(0.5)
        C = ((1/np.pi)**(0.5))*np.exp(-1*(X*10**-6)**2/4*kappa_p*t*10**-6)
        D = (-1*(X*10**-6)/(2*(kappa_p*t*10**-6)**0.5))*(erfc((X*10**-6)*(2*(kappa_p*t*10**-6)**0.5)))
        E = ((1/np.pi)**(0.5))*np.exp(-1*((X**10**-6)+b)**2/4*kappa_p*t*10**-6)
        F = (-1*((X*10**-6)**2+b**2)**0.5/(2*(kappa_p*t*10**-6)**0.5))*(erfc((((X*10**-6)**2+b**2)**0.5)/(2*(kappa_p*t*10**-6)**0.5)))
        return A*B*(C-D-E+F)

    #avoid mew and instead using q = Dn/t_c*A thus also encoding theta here


    temp = Temp(T) + 295 #to get in Kelvin
    soln.append(temp)
    Dissipated_E.append(e**2 * (m*g*h)*np.sin(theta)**2*D_frac(theta))
    max_T.append(np.amax(temp))
    
    """
    plt.figure()       
    CP = plt.contourf(X,T, temp, 100)
    cbar = plt.colorbar(CP, label ='Temperature/K')
    plt.title('Slab Temperature for Frictional Heating')
    plt.xlabel('distance within paper/10e-6m')
    plt.ylabel('time/10e-6s')
    plt.savefig('Friction_contour'+ str(np.degrees(theta)) +'.png')"""

plt.figure(figsize= [8,10])    
plt.subplot(4, 1, 1)
CP = plt.contourf(X,T, soln[1], 100)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Frictional Heating at ' + "{:.2f}".format(angle[1]) +' rad')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')

plt.subplot(4, 1, 2)
CP = plt.contourf(X,T, soln[4], 100)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Frictional Heating at ' + "{:.2f}".format(angle[4]) +' rad')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')

plt.subplot(4, 1, 3)
CP = plt.contourf(X,T, soln[6], 100)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Frictional Heating at ' + "{:.2f}".format(angle[6]) +' rad')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')

plt.subplot(4, 1, 4)
CP = plt.contourf(X,T, soln[8], 100)
cbar = plt.colorbar(CP, label ='Temperature/K')
plt.title('Slab Temperature for Frictional Heating at ' + "{:.2f}".format(angle[8]) +' rad')
plt.xlabel('distance within paper/10e-6m')
plt.ylabel('time/10e-6s')

plt.tight_layout()
plt.savefig('Friction_contour_subplots.png')

arrhenius = np.exp(-506/np.array(max_T))
    
angle_deg = np.degrees(angle)

plt.figure()
plt.plot(angle_deg,Dissipated_E)
plt.title('Dissipated Energy vs Angle of Incidence')
plt.xlabel('Angle of Incidence / degrees')
plt.ylabel('Dissipated Energy /  J')
plt.savefig('Friction_dissipated_E.png')

plt.figure()
plt.plot(angle_deg,arrhenius)
plt.title('Arrhenius Factor vs Angle of Incidence')
plt.xlabel('Angle of Incidence / degrees')
plt.ylabel('exp(- 506K / T)')
plt.savefig('Friction_dissipated_Arrhenius.png')