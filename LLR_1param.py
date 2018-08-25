#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 18:14:57 2018

@author: gabo
"""

import numpy as np
from math import factorial
from scipy.special import gammaln

def likelihood(Edad, rs, ts, e_m, R_s, R_f, tau=5730/np.log(2)):
    """Verosimilitud de una Edad de la muestra dadas las mediciones
    rs de la relación isotópica de la misma. Se asume la "fórmula de
    los alemanes" para la Edad y que toda la variabilidad es debida
    a la variabilidad en el conteo de 14C de la muestra.
    
    Esta función no es muy útil, porque los factoriales se vuelven tan
    grandes que no es posible convertirlos a float."""
    assert len(rs)==len(ts), "rs y ts deben tener igual longitud"
    R_m = R_s + (R_s - R_f) * np.exp(-Edad / tau)
    lambda_m = R_m * e_m
    L = 1
    for r, t in zip(rs, ts):
        mu = lambda_m * t
        k = np.floor(e_m * r) # Para evaluar el factorial, el arg debe ser int
        L = L * np.exp(-mu) * (mu)**(k) / factorial(k)
    return L

def loglikelihood(Edad, rs, ts, e_m, R_s, R_f, tau=5730/np.log(2)):
    """Loglikelihood de una Edad de la muestra dadas las mediciones
    rs de la relación isotópica de la misma. Se asume la "fórmula de
    los alemanes" para la Edad y que toda la variabilidad es debida
    a la variabilidad en el conteo de 14C de la muestra."""
    assert len(rs)==len(ts), "rs y ts deben tener igual longitud"
    R_m = R_s + (R_s - R_f) * np.exp(-Edad / tau)
    lambda_m = R_m * e_m
    # e_m * rs puede no ser entero, voy a usar la función gama en vez del
    # factorial.
    terminos = -lambda_m*ts + e_m*rs*np.log(lambda_m*ts) + gammaln(e_m*rs)
    return sum(terminos)
    

if __name__ == '__main__':
    
    charge_state = 3 
    q=charge_state*1.6e-10 # factor para pasar de nA a cargas
    tmedia = 5730. #periodo de semidesintegracion
    tau = tmedia/np.log(2.)
    
    #Muestra 
    c14m = np.array([159, 181, 133, 153, 169, 188, 205, 182, 169, 204, 215, 220, 242, 242, 216, 221, 235, 232, 308, 283, 257, 259], float)
    Iim = np.array([990, 923, 853, 862, 863, 	1030, 1030, 1035, 1070, 1130, 1160, 1140, 1180, 1180, 1200, 1200, 1250, 1200, 1230, 1250, 1320, 1220], float)
    Ifm = np.array([920, 854, 861, 864, 870, 1035, 1040, 1060, 1120, 1170, 1170, 1150, 1190, 1210, 1210, 1270, 1230, 1240, 1280, 1280, 1250, 1290], float)
    ts_m = 	300.00*np.ones(len(Iim))
    c12m = (Iim+Ifm)/2 * (ts_m/q)
    
    #Estándar
    c14s = np.array([160, 159, 179, 192, 203, 172, 202, 157],float)
    Iis = np.array([403, 406.5, 414, 416, 412, 415, 409, 411],float)
    Ifs	= np.array([406.5, 414, 415, 412, 414, 410, 410, 398],float)
    ts_s = 300.00*np.ones(len(Iis))
    c12s = (Iis+Ifs)/2 * (ts_s/q)
    r_s = c14s / c12s
    
    #Fondo 
    c14f = np.array([37, 25, 34, 53, 29, 32, 35, 25, 60], float)
    Iif = np.array([325.6, 366, 414, 585, 635, 673, 725, 749, 765], float)
    Iff = np.array([360, 416, 452, 628, 664, 695, 730, 760, 765], float)	#invente el ultimo
    ts_f = np.array([300, 300, 300, 300, 300, 300, 300, 300, 524], float)
    c12f = (Iif+Iff)/2 * (ts_f/q)
    
    
    rs_m = c14m / c12m
    e_m = np.average(c12m)
    R_s = np.average(c14s / c12s)
    R_f = np.average(c14f / c12f)

    L = lambda edad: loglikelihood(edad, rs_m, ts_m, e_m, R_s, R_f)
    # funciona!