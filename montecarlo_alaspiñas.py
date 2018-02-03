# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 18:00:46 2018

@author: Gabo
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import poisson, uniform

c14std, c12std_min, c12std_max = 178., 61621875000000., 61490625000000.
c14m, c12m_min, c12m_max = 212.4090909, 
c14f, c12f_min, c12f_max = 36.66666667, 96813333333333., 78416666666667.


N=1e4
rm = poisson.rvs(c14std, size=N)/uniform.rvs