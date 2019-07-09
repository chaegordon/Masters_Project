# -*- coding: utf-8 -*-
"""
Created on Sat May  4 13:25:37 2019

@author: chaeg
"""

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
    
df = pd.read_csv('data_1.csv', index_col = 0)

x_1 = df["x'"]
y_1 = df["y'"]


x_plot = np.linspace(0,1.5,100)

popt_1, pcov_1 = curve_fit(wiebull_1, x_1, y_1)

print("the fitting coefficients are")
print(popt_1)
print("the shapiro values are")
shapiro = stats.shapiro(y_1)
print(shapiro)
print("Anderson output is")
anderson = stats.anderson(y_1, dist = 'expon')
print(anderson)
print("The covariance matrix is")
print(pcov_1)


plt.figure()

plt.plot(x_plot, wiebull_1(x_plot, *popt_1), label = 'normal impact')
plt.title('Sphere on Sphere Collision')
plt.xlabel('Height/m')
plt.ylabel('Probability of Ignition')
plt.legend()
plt.savefig('prob_v_h_normal.png')