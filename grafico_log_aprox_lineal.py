# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 17:41:31 2018

@author: Samantha
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats


x = np.arange(0.001,10,0.01)
y = np.log(x)

inicio = 1.02
ancho_intervalo=0.5
int1, int2 = inicio, inicio+ancho_intervalo
ancho = 0.01
c1, c2 = np.mean([int1, int2])-ancho, np.mean([int1, int2])+ancho


x2 = np.arange(int1,int2,0.001)
y2 = np.log(x2)

x3 = np.arange(c1, c2, 0.001)
y3 = np.log(x3)

#x2 = x3
#y2 = y3

m, b = np.polyfit(x3,y3,1)

def recta(t):
    return m*t+b
    
r = recta(x)    
r2 = recta(x2)



fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y, 'b-', label= '$\log(x)$')
ax.plot(x,r, 'r-', label = 'Aprox. lineal')
ax.plot(x2, y2, 'k')
ax.plot(x2, r2, 'k')
#plt.plot(x, r, 'r')
#intervalos logaritmo
ax.vlines(x=int1, ymin=np.log(0.01),ymax=np.log(int1),colors='k', linestyles='dashed')
ax.vlines(x=int2, ymin=np.log(0.01),ymax=np.log(int2),colors='k', linestyles='dashed')
ax.hlines(y=np.log(int1), xmin=0,xmax=int1,colors='b', linestyles='dashed')
ax.hlines(y=np.log(int2), xmin=0,xmax=int2,colors='b', linestyles='dashed')
ax.vlines(x=0.011, ymin=np.log(int1), ymax=np.log(int2), colors='g', linestyles='solid')
#intervalos lineal
ax.vlines(x=int1, ymin=np.log(0.01),ymax=recta(int1),colors='k', linestyles='dashed')
ax.vlines(x=int2, ymin=np.log(0.01),ymax=recta(int2),colors='k', linestyles='dashed')
ax.hlines(y=recta(int1), xmin=0,xmax=int1,colors='r', linestyles='dashed')
ax.hlines(y=recta(int2), xmin=0,xmax=int2,colors='r', linestyles='dashed')
ax.vlines(x=0.011, ymin=recta(int1), ymax=recta(int2), colors='g', linestyles='solid')
ax.legend(loc=1, borderaxespad=0.3)
ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_xlim(xmin=0.01, xmax=0.5)
#ax.set_ylim(ymin=np.log(0.01), ymax=0)
ax.grid()

#plt.savefig('F:/Facultad/Estadística/loglineal.png', dpi=None, facecolor='w', edgecolor='w',         orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1, frameon=None)

#%%
def lineal(t,pendiente, ordenada):
    return pendiente*t+ordenada

ies = np.arange(0.01,1.5,0.01)
for i in ies:
    ancho_intervalo = 0.5
    int1_f, int2_f = i, i+ancho_intervalo
    ancho_f = 0.01
    c1_f, c2_f = np.mean([int1_f, int2_f])-ancho_f, np.mean([int1_f, int2_f])+ancho_f    
    x3_f = np.arange(c1_f, c2_f, 0.001)
    y3_f = np.log(x3_f)
    
    m_f, b_f = np.polyfit(x3_f,y3_f,1)
    intervalo_log = np.log(int2_f)- np.log(int1_f)
    intervalo_recta = lineal(int2_f, m_f, b_f)- lineal(int1_f, m_f, b_f)
    if intervalo_recta>intervalo_log:
        print(round(i,2), 'intervalo_recta>intervalo_log')
    elif intervalo_recta<intervalo_log:
        print(round(i,2), 'intervalo_log>intervalo_recta')
    
    
    
    
    
    