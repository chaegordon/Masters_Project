# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 16:12:07 2019

@author: chaeg
"""
v_2_0 = 0
t_1 = 0

import numpy as np

def v_1_2(theta, t):
    if t <= t_c:
      v_1_2 = v_1_0 - mew*(m_1/m_2)*v_2_0*(1-np.cos(omega*find_t2(theta)[0]))
      return v_1_2
    else:
      v_1_2 = v_1_0 - mew*(m_1/m_2)*v_2_0*(1-np.cos(((omega/e)*find_t2(theta)[0])+(np.pi/2)*(1 - 1/e))
      return v_1_2
        