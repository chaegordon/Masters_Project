# -*- coding: utf-8 -*-
"""
Created on Mon May 13 00:02:00 2019

@author: chaeg
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from scipy import stats

df = pd.read_csv('data_V.csv', index_col = 0)
df = df.loc[(df["v1"]!=0)]
df = df.loc[(df["v1"]!=1)]
df = df.loc[(df["P"]!=0)]
df = df.loc[(df["P"]!=1)]

x_1 = df["v1"]
x_err = 0.0011*np.ones(x_1.shape)
y_1 = np.array(df["P"])
y_err1 = (y_1*(1-y_1)/len(y_1))**(0.5)

# Linear regression using stats.linregress
slope, intercept, r_value, p_value, std_err = stats.linregress(x_1, y_1)
print('Linear regression using stats.linregress Energy')
print("slope: %f    intercept: %f" % (slope, intercept))
print("r-squared: %f" % r_value**2)
print(std_err)


plt.figure()
plt.scatter(x_1,y_1, label = 'Burn-Velocity Correlation')
plt.plot(np.unique(x_1), np.poly1d(np.polyfit(x_1, y_1, 1))(np.unique(x_1)))
plt.errorbar(x_1,y_1, yerr=y_err1, xerr=x_err, fmt = "none")
plt.title('Burn Probability Against Velocity')
plt.xlabel(r'Vertical Velocity of Impact/ms^-1')
plt.ylabel(r'Burn Probability')
plt.legend(loc='best')
plt.savefig('Appendix_Burn_V_intercept.pdf',Bbox = 'tight')