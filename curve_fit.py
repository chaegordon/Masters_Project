# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 20:28:14 2019

@author: chaeg
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from scipy import stats


def func(x,a,b):
   return a*x - a*np.log(b)

"""by taking logs of the functions the data can be linearised as above however,
because of 0's in the logs the points on the plateaus of the distribution are
not defined and so this method is not suitable to my project due to the small
number of remaining data points"""

#if curve fit error less than 1/3 of statistical error then no biggie

#Need to use  the coefficients to calculate the errors in the data

def wiebull(x,a,b):
    return 1 - np.exp(-1*(x/b)**a)

#need to renormalise for the different maximum (factor of 0.8), could have included another parameter in fit
def wiebull_1(x,a,b):
    return 0.8 - 0.8*np.exp(-1*(x/b)**a)

def wiebull_2(x,a,b):
    return 0.28 - 0.28*np.exp(-1*(x/b)**a)
    
def h_50(a,b):
    return (b*(np.log(2))**(1/a))**-1

def h_50_1(a,b):
    return (b*(np.log(0.8/0.3))**(1/a))**-1

def h_50_2(a,b):
    return 0
    
df = pd.read_csv('data_1.csv', index_col = 0)

x_1 = df["x'"]
y_1 = df["y'"]

x_2 = df["x''"]
y_2 = df["y''"]

x_3 = df["x'''"]
y_3 = df["y'''"]

x_4 = df["x4"]
y_4 = df["y4"]

x_5 = df["x5"]
y_5 = df["y5"]

x_6 = df["x6"]
y_6 = df["y6"]

x_plot = np.linspace(0,1.5,100)

popt_1, pcov_1 = curve_fit(wiebull_1, x_1, y_1)

popt_2, pcov_2 = curve_fit(wiebull, x_2, y_2)

popt_3, pcov_3 = curve_fit(wiebull, x_3, y_3)

popt_4, pcov_4 = curve_fit(wiebull_2, x_4, y_4)

popt_5, pcov_5 = curve_fit(wiebull, x_5, y_5)

popt_6, pcov_6 = curve_fit(wiebull, x_6, y_6)

h_theta = []
h_theta.append(h_50_1(*popt_1))
h_theta.append(h_50(*popt_2))
h_theta.append(h_50(*popt_3))
h_theta.append(h_50_2(*popt_4))
h_theta.append(h_50(*popt_5))
h_theta.append(h_50(*popt_6))

h_df = pd.DataFrame(h_theta,columns = ['1/h_50'])

""" rename the angle so it is consistent with the literature definition (90-angle)"""

theta_df = pd.DataFrame([90,29,45,17,60,70], columns =['Angle of Incidence/ degrees'])

h_t_df = theta_df.join(h_df)



plt.figure()
ax = h_t_df.plot.scatter(x='Angle of Incidence/ degrees', y='1/h_50')
fig = ax.get_figure()
fig.suptitle('1/h_50 versus Angle of Incidence', fontsize=16)
fig.savefig('h_50_versus_theta.pdf')

#print(popt_1)
#CP = plt.scatter(x,y)
"""Also renormalised here or can use wiebull_1 in the plot same difference"""

print(popt_1)

shapiro = stats.shapiro(y_1)
print(shapiro)

plt.figure()

plt.plot(x_plot, wiebull_1(x_plot, *popt_1), label = 'normal impact')

plt.plot(x_plot, wiebull(x_plot, *popt_2), label = '29 degrees impact')

plt.plot(x_plot, wiebull(x_plot, *popt_3), label = '45 degrees impact')

plt.plot(x_plot, wiebull_2(x_plot, *popt_4), label = '17 degrees impact')

plt.plot(x_plot, wiebull(x_plot, *popt_5), label = '60 degrees impact')

plt.plot(x_plot, wiebull(x_plot, *popt_6), label = '70 degrees impact')

plt.title('Drop Test Data')
plt.xlabel('Height/m')
plt.ylabel('Probability of Ignition')
plt.legend()
plt.savefig('prob_v_h.pdf')