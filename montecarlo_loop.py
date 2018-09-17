#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 18:50:50 2018

@author: gabo
"""

import numpy as np
from scipy.stats import poisson, uniform, beta, gamma
import matplotlib.pyplot as plt
import time
#%%
#Levanto los datos
#corrientes en nA y tiempos en segundos

#STD 
Iistd = np.array([403, 406.5, 414, 416, 412, 415, 409, 411],float)
Ifstd	= np.array([406.5, 414, 415, 412, 414, 410, 410, 398],float)
c14std = np.array([160, 159, 179, 192, 203, 172, 202, 157],float)
dtstd = 300.00*np.ones(len(Iistd))

#Muestra 
Iim = np.array([990, 923, 853, 862, 863, 	1030, 1030, 1035, 1070, 1130, 1160, 1140, 1180, 1180, 1200, 1200, 1250, 1200, 1230, 1250, 1320, 1220], float)
Ifm = np.array([920, 854, 861, 864, 870, 1035, 1040, 1060, 1120, 1170, 1170, 1150, 1190, 1210, 1210, 1270, 1230, 1240, 1280, 1280, 1250, 1290], float)
c14m = np.array([	159, 181, 133, 153, 169, 188, 205, 182, 169, 204, 215, 220, 242, 242, 216, 221, 235, 232, 308, 283, 257, 259], float)
dtm = 	300.00*np.ones(len(Iim))

#Fondo 
Iif = np.array([325.6, 366, 414, 585, 635, 673, 725, 749, 765], float)
Iff = np.array([360, 416, 452, 628, 664, 695, 730, 760, 765], float)	#invente el ultimo
c14f = np.array([	37, 25, 34, 53, 29, 32, 35, 25, 60], float)
dtf = np.array([300, 300, 300, 300, 300, 300, 300, 300, 524], float)

charge_state = 3 
q=charge_state*1.6e-10 # factor para pasar de nA a cargas
tmedia = 5730. #periodo de semidesintegracion
tau = tmedia/np.log(2.)

# Conteos de 12C iniciales y finales para el estandar, muestra y fondo.
c12istd, c12fstd = Iistd*dtstd/q, Ifstd*dtstd/q 
c12im, c12fm = Iim*dtm/q, Ifm*dtm/q 
c12if, c12ff = Iif*dtf/q, Iff*dtf/q 

#%%

#Tomo los promedios de los conteos de 14C
c14std_p, c14m_p, c14f_p = np.mean(c14std), np.mean(c14m), np.mean(c14f)
print("--- %s cuentas de 14C durante la medición de la muestra" % c14m_p)
print("--- %s cuentas de 14C durante la medición del fondo" % c14f_p)
print("--- %s cuentas de 14C durante la medición del estandar" % c14std_p)

c12im_p, c12fm_p = np.mean(c12im), np.mean(c12fm)
c12m=(c12im_p+c12fm_p)/2
print("--- %s cuentas de 12C durante la medición de la muestra" % c12m)

c12if_p, c12ff_p = np.mean(c12if), np.mean(c12ff)
c12f=(c12if_p+c12ff_p)/2
print("--- %s cuentas de 12C durante la medición del fondo" % c12f)

c12istd_p, c12fstd_p = np.mean(c12istd), np.mean(c12fstd)
c12std=(c12istd_p+c12fstd_p)/2
print("--- %s cuentas de 12C durante la medición del estandar" % c12std)

#%%
# Ahora reemplazamos la suma de todos los c14m por otros valores posibles,
# no muy lejanos al original (4673)
n_edades = 100
cant = int(1e6) # cantidad de nros random en cada simu

#sumas_c14m = np.geomspace(1700, 6000, n_edades)
sumas_c14m = np.linspace(1700, 6000, n_edades)
#%%
start_time = time.time()
intervalos_edad = np.zeros((n_edades, 2), dtype=int)
edades_medias, errores_relativos = np.zeros((2, n_edades))

for i, suma_c14m in enumerate(sumas_c14m):
    c14m_sim = poisson.rvs(suma_c14m, size=cant) / np.array(len(c14m), float)
    c14f_sim = poisson.rvs(np.sum(c14f), size=cant) / np.array(len(c14f), float)
    c14std_sim = poisson.rvs(np.sum(c14std), size=cant) / np.array(len(c14std), float)
    
    Rm = c14m_sim/uniform.rvs(loc=c12im_p, scale=c12fm_p-c12im_p, size=cant)
    Rf = c14f_sim/uniform.rvs(loc=c12if_p, scale=c12ff_p-c12if_p, size=cant)
#     Para el estándar, la corriente final es más chica que la inicial
    # Por lo tanto los extremos de la uniforme son al revés que en los otros casos
    Rstd = c14std_sim/uniform.rvs(loc=c12fstd_p, scale=c12istd_p-c12fstd_p, size=cant)
    
    # Testing con 12C fijos
#    Rm = c14m_sim / ((c12im_p + c12fm_p) / 2)
#    Rf = c14f_sim/ ((c12if_p + c12ff_p) / 2)
#    Rstd = c14std_sim / ((c12fstd_p + c12istd_p) / 2)
    
    Rtot = (Rm-Rf)/(Rstd-Rf)
    
    #Elimino los infinitos
    mask = (Rtot != np.float('+inf'))
    Rtot = Rtot[mask]
    # Elimino los posibles casos en que esta variable tome valores negativos
    mask = (Rtot>=0.00)
    Rtot = Rtot[mask]
    edad = -tau*np.log(Rtot)
#    import pdb; pdb.set_trace()
    l, u = np.percentile(edad, 16), np.percentile(edad, 84) #lower, upper
    intervalos_edad[i,:] = np.array([l, u], dtype=int)
    edades_medias[i] = (l + u) / 2
    errores_relativos[i] = ((u - l) / 2) / edades_medias[i]
#    import pdb; pdb.set_trace()
print("--- %s seconds ---" % (time.time() - start_time))
#%%
with plt.style.context(('seaborn')):
    fig, ax = plt.subplots()
ax.plot(edades_medias, errores_relativos * 100, '--.k')
ax.set_xlabel('Edad (años)', fontsize=18)
ax.set_ylabel('Error relativo (%)', fontsize=18)
ax.tick_params(labelsize=16)
#
#ax.axvline(8482, color='r', label='Valor CEMA')
#ax.legend(fontsize=18)

fig.tight_layout()

#%%
with plt.style.context(('seaborn')):
    fig, ax = plt.subplots()
for i in range(n_edades):
    ax.plot([intervalos_edad[i,0], intervalos_edad[i,1]], [sumas_c14m[i]/22] * 2, '-c')
ax.set_xlabel('Edad (años)', fontsize=18)
ax.set_ylabel(r'$^{14}C$ de la muestra (c)', fontsize=18)
ax.tick_params(labelsize=16)
#
#ax.axvline(8482, color='r', label='Valor CEMA')
#ax.legend(fontsize=18)

fig.tight_layout()
#%%
with plt.style.context(('seaborn')):
    fig, ax = plt.subplots()
ax.plot(edades_medias, edades_medias, 'k-')
ax.plot(edades_medias, intervalos_edad[:,0], '-c')
ax.plot(edades_medias, intervalos_edad[:,1], '-c')
ax.set_xlabel('Edad (años)', fontsize=18)
ax.set_ylabel('Intervalo 68% C.L. (años)', fontsize=18)
ax.tick_params(labelsize=16)
#
#ax.axvline(8482, color='r', label='Valor CEMA')
#ax.legend(fontsize=18)

fig.tight_layout()