# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 17:50:25 2018

@author: Samantha
"""

import numpy as np
from sympy import *
from sympy import init_printing
init_printing() 

#Variables
c14m, c14f,c14s, i_m, i_f, i_s, tiempo, tau, carga = symbols('C_{14M} C_{14F} C_{14S} I_M I_F I_S t T q')
#Errores
dcm, dcf, dcs, dim, dif, dis = symbols('e_{CM} e_{CF} e_{CS} e_{IM} e_{IF} e_{IS}')


#Definimos las funciones de interes
def c12(i):
    return (i*tiempo)/carga
c12m, c12f, c12s = c12(i_m), c12(i_f), c12(i_s)
R_m = c14m/c12m
R_f = c14f/c12f
R_s = c14s/c12s
edad = -tau * log((R_m - R_f) / (R_s - R_f))

edad
#%%

#Derivo
dedc14m = diff(edad, c14m)
dedc14f = diff(edad, c14f)
dedc14s = diff(edad, c14s)
dedim = diff(edad, i_m)
dedif = diff(edad, i_f)
dedis = diff(edad, i_s)

error_edad = sqrt((dedc14m**2)*(dcm**2) + (dedc14f**2)*(dcf**2) + (dedc14s**2)*(dcs**2)  + (dedim**2)*(dim**2) + (dedif**2)*(dif**2) + (dedis**2)*(dis**2))

error_edad

#%%
#Reemplazo con los valores

#Levanto los datos
#corrientes en nA y tiempos en segundos

#STD 
Iistd_medicion = np.array([403, 406.5, 414, 416, 412, 415, 409, 411],float)
Ifstd_medicion	= np.array([406.5, 414, 415, 412, 414, 410, 410, 398],float)
c14std_medicion = np.array([160, 159, 179, 192, 203, 172, 202, 157],float)
dtstd_medicion = 300.00*np.ones(len(Iistd_medicion))

#Muestra 
Iim_medicion = np.array([990, 923, 853, 862, 863, 	1030, 1030, 1035, 1070, 1130, 1160, 1140, 1180, 1180, 1200, 1200, 1250, 1200, 1230, 1250, 1320, 1220], float)
Ifm_medicion = np.array([920, 854, 861, 864, 870, 1035, 1040, 1060, 1120, 1170, 1170, 1150, 1190, 1210, 1210, 1270, 1230, 1240, 1280, 1280, 1250, 1290], float)
c14m_medicion = np.array([	159, 181, 133, 153, 169, 188, 205, 182, 169, 204, 215, 220, 242, 242, 216, 221, 235, 232, 308, 283, 257, 259], float)
dtm_medicion = 	300.00*np.ones(len(Iim_medicion))

#Fondo 
Iif_medicion = np.array([325.6, 366, 414, 585, 635, 673, 725, 749, 765], float)
Iff_medicion = np.array([360, 416, 452, 628, 664, 695, 730, 760, 765], float)	#invente el ultimo
c14f_medicion = np.array([	37, 25, 34, 53, 29, 32, 35, 25, 60], float)
dtf_medicion = np.array([300, 300, 300, 300, 300, 300, 300, 300, 524], float)

charge_state = 3 
q=charge_state*1.6e-10 # factor para pasar de nA a cargas
tmedia = 5730. #periodo de semidesintegracion
tau_valor = tmedia/np.log(2.)

# Lo que sigue es válido sólo cuando las mediciones son todas de igual duración.
# Aunque no es necesario tener el mismo número de mediciones para todos los casos.

#Tomo los promedios de los conteos de 14C (esto en realidad no esta bien)
c14std_p, c14m_p, c14f_p = np.mean(c14std_medicion), np.mean(c14m_medicion), np.mean(c14f_medicion)

i_m_p = (np.max(Ifm_medicion)+np.min(Iim_medicion))/2
i_f_p = (np.max(Iff_medicion)+np.min(Iif_medicion))/2
i_s_p = (np.max(Ifstd_medicion)+np.min(Iistd_medicion))/2

err_c14std_p, err_c14m_p, err_c14f_p = np.sqrt(c14std_p), np.sqrt(c14m_p), np.sqrt(c14f_p)

err_i_m_p = (np.max(Ifm_medicion)-np.min(Iim_medicion))/2
err_i_f_p = (np.max(Iff_medicion)-np.min(Iif_medicion))/2
err_i_s_p = (np.max(Ifstd_medicion)-np.min(Iistd_medicion))/2

#%%

#Numeros
edad_valor = edad.subs({c14m:c14m_p, c14f:c14f_p,c14s:c14std_p, i_m:i_m_p, i_f:i_f_p, i_s:i_s_p, tiempo:300, tau:tau_valor, carga:q})

edad_valor

#%%

error_edad_valor = error_edad.subs({c14m:c14m_p, c14f:c14f_p,c14s:c14std_p, i_m:i_m_p, i_f:i_f_p, i_s:i_s_p, tiempo:300, tau:tau_valor, carga:q, dcm:err_c14m_p, dcf:err_c14f_p, dcs:err_c14std_p, dim:err_i_m_p, dif:err_i_f_p, dis:err_i_s_p})

error_edad_valor

#%%

#Si solo consideramos que varia el carbono 14 de la muestra

error_edad_c14m = sqrt((dedc14m**2)*(dcm**2))

error_edad_c14m_valor = error_edad_c14m.subs({c14m:c14m_p, c14f:c14f_p,c14s:c14std_p, i_m:i_m_p, i_f:i_f_p, i_s:i_s_p, tiempo:300, tau:tau_valor, carga:q, dcm:err_c14m_p, dcf:err_c14f_p, dcs:err_c14std_p, dim:err_i_m_p, dif:err_i_f_p, dis:err_i_s_p})

error_edad_c14m_valor



