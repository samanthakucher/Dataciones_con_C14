# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:53:45 2018

@author: Samantha
"""

import numpy as np
from scipy.stats import poisson, uniform
import matplotlib.pyplot as plt
import time

#Levanto los datos

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

q = 1e-9 #estado de carga
tmedia = 5730. #semi periodo de desintegracion
tau = tmedia/np.log(2.)

c12istd, c12fstd = Iistd*dtstd/q, Ifstd*dtstd/q
c12im, c12fm = Iim*dtm/q, Ifm*dtm/q
c12if, c12ff = Iif*dtf/q, Iff*dtf/q

#%%

#Tomo los promedios

c14std_p, c14m_p, c14f_p = np.mean(c14std), np.mean(c14m), np.mean(c14f)

c12istd_p, c12fstd_p = np.mean(c12istd), np.mean(c12fstd)
c12im_p, c12fm_p = np.mean(c12im), np.mean(c12fm)
c12if_p, c12ff_p = np.mean(c12if), np.mean(c12ff)
#%%

#Cantidad de valores que voy a generar en cada Montecarlo
cant = int(1e6)

start_time = time.time()

Rstd = poisson.rvs(c14std_p, size=cant)/uniform.rvs(loc=c12istd_p, scale=c12fstd_p, size=cant)
Rm = poisson.rvs(c14m_p, size=cant)/uniform.rvs(loc=c12im_p, scale=c12fm_p, size=cant)
Rf = poisson.rvs(c14f_p, size=cant)/uniform.rvs(loc=c12if_p, scale=c12ff_p, size=cant)

Rtot = (Rm-Rf)/(Rstd-Rf)


#Elimino los infinitos
mask = (Rtot != np.float('+inf'))
Rtot = Rtot[mask]

# Verifico que para todos los elementos sea cierto que no son +inf
np.all(Rtot != np.float('+inf'))

print("--- %s seconds ---" % (time.time() - start_time))
#%%

Rmin, Rmax = 0, 1.25

puntos_a = np.linspace(Rmin, Rmax,300) #puntos para graficar f(x)
bines = np.linspace(Rmin, Rmax,100)
numero, bins = np.histogram(Rtot, bins = bines) #numero=numero de entradas por bin
error = np.sqrt(numero) / (np.diff(bins)* np.sum(numero)) #error poissoniano
numero = numero / (np.diff(bins) * np.sum(numero)) #Normalizo a 1 (divido por el área ocupada por el histograma)

fig = plt.figure(figsize=(10,6))
plt.bar(bins[:-1], numero, width = np.diff(bins), yerr = error, ecolor="b", color='c', alpha=0.7)
plt.legend(loc=1, borderaxespad=0.)
plt.xlim([Rmin, Rmax])
plt.xlabel('R')
#plt.ylabel('f(x)')
plt.title('Histograma')
plt.grid()
plt.show()

#%%

edad = -tau*np.log(Rtot)
#Acá marca error porque divide por cero (-log(1))
