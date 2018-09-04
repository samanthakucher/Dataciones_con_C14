# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 17:50:25 2018

@author: Samantha
"""

import numpy as np
from sympy import *
from sympy import init_printing
init_printing() 

charge_state = 3 
q=charge_state*1.6e-10 # factor para pasar de nA a cargas
tmedia = 5730. #periodo de semidesintegracion
#tau = tmedia/np.log(2.)
      
c14m, c14f,c14s, i_m, i_f, i_s, tiempo, tau, carga = symbols('C^{14}_M C^{14}_F C^{14}_S I_M I_F I_S t T q')

#dc14m, dc14f,dc14s, diim, difm, diif, diff, diis, difs = symbols('DC^{14}_M DC^{14}_F DC^{14}_S DI^i_M DI^f_M DI^i_F DI^f_F DI^i_S DI^f_S')

def c12(i):
    return (i*tiempo)/carga

c12m, c12f, c12s = c12(i_m), c12(i_f), c12(i_s)

R_m = c14m/c12m
R_f = c14f/c12f
R_s = c14s/c12s
edad = -tau * log((R_m - R_f) / (R_s - R_f))

edad
#%%
dedc14m = diff(edad, c14m)
dedc14f = diff(edad, c14f)
dedc14s = diff(edad, c14s)
dedim = diff(edad, i_m)
dedif = diff(edad, i_f)
dedis = diff(edad, i_s)



