# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 00:29:22 2018

@author: Samantha
"""

import numpy as np
from scipy.stats import uniform
import matplotlib.pyplot as plt

M=int(1e6)
x = uniform.rvs(size=M) #numeros aleatorios con distribucion uniforme
y = x**2

#%%

bines = np.linspace(0,1.1,300)
numero1, bins1 = np.histogram(x, bins = bines) #numero=numero de entradas por bin
numero2, bins2 = np.histogram(y, bins = bines)
#error = np.sqrt(numero) / (np.diff(bins)* np.sum(numero)) #error poissoniano
numero1 = numero1 / (np.diff(bins1) * np.sum(numero1)) #Normalizo a 1 (divido por el área ocupada por el histograma)
numero2 = numero2 / (np.diff(bins2) * np.sum(numero2))

fig = plt.figure(figsize=(10,2))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
ax1.bar(bins1[:-1], numero1, width = np.diff(bins1), color='c', alpha=0.5)
ax1.set_xlim([0,1])
ax1.set_xlabel('x', fontsize=12)
ax1.set_ylabel('f(x)', fontsize=12)
#plt.title('$x$')
ax1.grid()
ax2.bar(bins2[:-1], numero2, width = np.diff(bins2), color='c', alpha=0.5)
ax2.set_xlim([0,1])
ax2.set_ylim([0,7])
ax2.set_xlabel('y', fontsize=12)
ax2.set_ylabel('f(y)', fontsize=12)
#plt.title('$x$')
ax2.grid()
plt.show()

plt.savefig('F:/x_xcuadrado.jpg', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', frameon=None)