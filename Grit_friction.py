import numpy as np
import matplotlib.pyplot as plt
from math import factorial

"""estimate rho, estimated p as 0.5"""

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
Area = np.pi*(b/2)**2
tau = 400*10**3 # using burst strength of paper
x_0 = 6.75*10**-5 #max comp. steel on steel sphere


l = 100 #micrometers is the paper width
N = np.linspace(0,4,5,dtype = float)
n = np.linspace(0,99,100,dtype = int)
x=np.linspace(0*l,l,10**2)
time=np.logspace(-2,-1,10**3, dtype = float) #in microseconds
X,T = np.meshgrid(x,time)

R = 0.05
M = 0.5

#NEED t-C calculated each time
t_c = 5*10**-6

#material stiffness constants
E_m = 200*10**9


angle = np.linspace(0.1,1.56,50)

#get w by fitting data
w=125600
e= 0.2

soln_grit = []
soln_sub = []

for theta in angle:
    #need to make this E-eff
    t_c = (1+e)*2.87*(((M/2)**2)/float(R*E_m**2*V*np.cos(theta))**0.2)
    A_tot = np.sin(theta)*2*np.pi*R*V/w + 2*t_c*np.cos(theta)*V**(1.5)*(2*R*np.sin(theta)/w)**0.5
    L=t_c*V*np.cos(theta)
    p=0.5
    rho = 0.1*(10*10**6)**2 #number per m^2, was 220*10**9, new value taken from SEM from paper
    
    def P_gg(k,j):
        A=((rho*A_tot)**k)/factorial(k)
        B= np.exp(-1*rho*A_tot)
        C=factorial(k)//(factorial(j)*factorial(k-j))
        D= (p**j)*(1-p)**(k-j)
        E= (1-((A_tot/float(4*R*L*(k-j+1)))*(1-(1-(4*R*L/float(A_tot))**(k-j+1))))**j)
        return A*B*C*D*E
    
    def P_tot_gg(n):
        counter = []
        for k in range(1,n):
            for j in range(1,k):
                counter += P_gg(k,j)
        return sum(counter)
    
    P_tot_grit = P_tot_gg(20)           
        
    P_tot_substrate = 1 - np.exp(-1*rho*A_tot)
    soln_grit.append(P_tot_grit)
    soln_sub.append(P_tot_substrate)
    
    
angle_deg = np.degrees(angle)
plt.figure()
plt.scatter(angle_deg,soln_grit)
plt.title('Probability vs Angle of Incidence for Grit on Grit Friction')
plt.xlabel('Angle of Incidence/degrees')
plt.ylabel('Probability')
plt.savefig('prob_grit.png',bbox_inches='tight')

plt.figure()
plt.scatter(angle_deg,soln_sub)
plt.title('Probability vs Angle of Incidence for Grit on Substrate Friction')
plt.xlabel('Angle of Incidence/degrees')
plt.ylabel('Probability')
plt.savefig('prob_sub.png',bbox_inches='tight')    
#shows that grit could be the cause of the ignition however, there must be some
#more dynamical details?
    