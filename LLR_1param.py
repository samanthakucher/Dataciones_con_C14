#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 18:14:57 2018

@author: gabo
"""

import numpy as np
from scipy.special import gammaln
import matplotlib.pyplot as plt

def Edad_func(R_m, R_s, R_f, tau=5730/np.log(2)):
    return -tau * np.log((R_m - R_f) / (R_s - R_f))

def R_m_func(Edad, R_f, R_s, tau=5730/np.log(2)):
    return R_f + (R_s - R_f) * np.exp(-Edad / tau)

def loglikelihood(Edad, rs, ts, e_m, R_s, R_f, tau=5730/np.log(2)):
    """Loglikelihood de una Edad de la muestra dadas las mediciones
    rs de la relación isotópica de la misma. Se asume la "fórmula de
    los alemanes" para la Edad y que toda la variabilidad es debida
    a la variabilidad en el conteo de 14C de la muestra.
    Si rs, ts son arrays 2d, para poder evaluar la función en múltiples puntos,
    la primera dimensión debe indexar el vector en el que se quiere evaluar y
    la segunda dimensión debe indexar las componentes del vector. Edad también
    puede ser un array en vez de un float (pero solo 1d)."""
    assert rs.shape == ts.shape, "rs y ts deben tener igual shape"
    assert (rs.ndim == 1) or (rs.ndim == 2), "rs, ts deben ser 1d, o 2d"
    
    R_m = R_m_func(Edad, R_f, R_s)
    lambda_m = R_m * e_m
    
    if isinstance(Edad, np.ndarray): # Para broadcasting correcto
        if rs.ndim == 1:
            lambda_m = lambda_m[:,np.newaxis]
        elif rs.ndim == 2:
            lambda_m = lambda_m[:,np.newaxis, np.newaxis]
            
    # e_m * rs puede no ser entero, voy a usar la función gama en vez del
    # factorial. Recordar que si n es entero positivo, n! == gamma(n+1)
    # y que gamma para los reales negativos es otra cosa totalmente distinta
    terminos = ( -lambda_m * ts + e_m * rs * ts * np.log(lambda_m * ts)
                - gammaln(e_m * rs * ts + 1) )
    return np.sum(terminos, axis=-1)
#    return terminos[:, 0] # Prueba

    

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


    LL_fijtot = lambda edad: loglikelihood(edad, rs_m, ts_m, e_m, R_s, R_f)
    # funciona!
#%%
    # Veamos que se obtiene de reemplazar estos datos en la fórmula para la edad
    edad_esperada = Edad_func(rs_m, R_s, R_f)
    print('Edades que se obtienen con cada valor de r_m, dejando fijos R_s y R_f:')
    print(('{:.0f} '*len(rs_m)).format(*edad_esperada))
    print('Promedio y desviación estándar: {:.0f} +/- {:.0f} años'.format(np.mean(edad_esperada),
          np.std(edad_esperada)))

    # Es decir que con estos datos deberíamos ver que la verosimilitud
    # se maximiza alrededor de Edad = 9000
#%%
    edades1 = np.linspace(0, int(1e4), 1000)
    lls1 = LL_fijtot(edades1)
    edades2 = np.linspace(0, int(1e5), 10000)
    lls2 = LL_fijtot(edades2)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(edades1, lls1, '-', color='dodgerblue')
    ax2.plot(edades2, lls2, '-', color='dodgerblue')
    fig.tight_layout()
    
    # Vemos que el loglikelihood tiene un máximo aprox en 9000 años.
    # Hacemos la cuenta de la edad de Máxima Verosimilitud de manera analítica:
    
    r_bar = np.average(rs_m)
    edad_MV = Edad_func(r_bar, R_s, R_f)
    print('La edad de Máxima Verosimilitud es {:.0f} años'.format(edad_MV))
    # La edad de Máxima Verosimilitud es 8797 años
#%%
    # Hago una grilla cerquita de ese valor p/ determinar el máximo
    # con mayor precisión.
    es = np.linspace(8000, 10000, int(1e6))
    precision_horizontal = np.diff(es)[0] # error debido al mesh despreciable
    lls = LL_fijtot(es)
    llmax = np.max(lls)
    ind_intervalo = np.where(lls >= (llmax - 1/2))[0]
    ind_L, ind_U = min(ind_intervalo), max(ind_intervalo)
    e_L, e_U = es[ind_L], es[ind_U]
    argllmax = es[np.argmax(lls)]
    
    with plt.style.context(('seaborn')):
        fig, ax = plt.subplots()
    ax.plot(es, lls, '-', color='dodgerblue')
    ax.plot(argllmax, llmax, 'ro')
    ax.axvline(e_L, color='deeppink')
    ax.axvline(e_U, color='deeppink')
    print('Gráficamente, la edad de Máxima Verosimilitud es {:.0f} años.'.format(argllmax))
    print('El intervalo LLR con esta verosimilitud es [{:.0f}, {:.0f}] años.'.format(e_L, e_U))
    print('Error relativo del {:.2f}%.'.format((e_U - e_L) / argllmax * 100))
    
#%%    
    # Por otro lado, veamos si pasa algo razonable al considerar la verosim
    # como función de los datos, con parámetro fijo.
    edad = edad_MV
    # El valor para el que se maximiza la prob. debería ser alrededor de 
    rmax = R_m_func(edad, R_f, R_s) # -> 3.040115467990185e-13
    # Sup que todas las 22 mediciones dan exactamente lo mismo, y variamos
    # el valor que pueden adoptar:
    
    def variar_valor_rs(npuntos, rmin, rmax, tiempo):
        ts = np.ones((npuntos, 22)) * tiempo # tiempo de medición
        rs = np.zeros((npuntos, 22)) 
        for i in range(22):
            rs[:,i] = np.linspace(rmin, rmax, npuntos)
        return ts, rs
    
    LL_fij3 = lambda edad, rs_m, ts_m: loglikelihood(edad, rs_m, ts_m, e_m, R_s, R_f)
    
    ts, rs = variar_valor_rs(int(1e5), 0, 1e-12, 300) # meds de 300 segundos
    lls = LL_fij3(edad, rs, ts)
    fig, ax = plt.subplots()
    ax.plot(rs[:,0], lls, '-', color='deeppink')
    
    # Efectivamente, vemos un máximo en el lugar predicho.
    # Ya que estamos, grafiquemos la probabilidad también:
    
    ts, rs = variar_valor_rs(int(1e5), 1e-13, 5e-13, 300)
    lls = LL_fij3(edad, rs, ts)
    ls = np.exp(lls)
    fig, ax = plt.subplots()
    ax.plot(rs[:,0], ls, '-', color='deeppink')
    
    # Guau, la distribución es básicamente una gaussiana muy, muy angosta,
    # de sigma ~ 1e-15.

    