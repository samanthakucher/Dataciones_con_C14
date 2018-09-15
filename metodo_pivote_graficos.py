# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

import numpy as np
from scipy.stats import binom, poisson, norm
from matplotlib import pyplot as plt

def find_index_of_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def mover_mu(valor): #en valor tiene que ir el mu que elegimos para graficar
    integral = np.zeros(shape=int(2*valor))
    for i in range(0,int(2*valor)):
        integral[i] = poisson.sf(k=i, mu=valor, loc=0)

    mu_max = find_index_of_nearest(integral, 0.16)
    mu_min = find_index_of_nearest(integral, 0.84)

    return mu_min, mu_max
    
valor_elegido = 20
mu_max, mu_min = mover_mu(valor_elegido)

x = np.arange(0,int(2*valor_elegido))
p1, p2, p3 = poisson.pmf(x,valor_elegido), poisson.pmf(x,mu_min), poisson.pmf(x,mu_max)

fig = plt.figure(figsize=(7,7))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

fig.text(0.04, 0.5, 'P(k)', va='center', rotation='vertical')

ax1.plot(x, p1, 'b-', label = 'Poisson($\mu$=20)')
ax1.axvline(x=valor_elegido, color='r')
ax1.set_ylim(ymax=0.12)
ax1.set_yticks([0,0.05,0.1])
#ax1.legend(loc=0, borderaxespad=0.5)
#ax1.set_xlabel('k')
#ax1.set_ylabel('P(k)')
#ax1.set_title('Grafico 1')
ax1.grid()

ax2.plot(x, p2, 'b-', label = 'Poisson($\mu$=24)')
ax2.set_ylim(ymax=0.12)
ax2.set_yticks([0,0.05,0.1])
ax2.axvline(x=valor_elegido, color='r')
ax2.fill_between(x = x[:valor_elegido+1], y1 = p2[:valor_elegido+1], color='c', label='$16\%$')
#ax2.legend(loc=0, borderaxespad=0.5)
#ax2.set_xlabel('k')
#ax2.set_ylabel('P(k)')
#ax3.set_title('Grafico 1')
ax2.grid()

ax3.plot(x, p3, 'b-',label = 'Poisson($\mu$=15)')
ax3.set_ylim(ymax=0.12)
ax2.set_yticks([0,0.05,0.1])
ax3.axvline(x=valor_elegido, color='r')
ax3.fill_between(x = x[valor_elegido:], y1 = p3[valor_elegido:], color='c', label='$16\%$')
ax3.set_xlabel('k')
#ax3.legend(loc=0, borderaxespad=0.5)
#ax3.set_ylabel('P(k)')
#ax3.set_title('Grafico 1')
ax3.grid()

#plt.savefig('F:/Facultad/Estadística/pivote_sinlabels.png', dpi=None, facecolor='w', edgecolor='w',         orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1, frameon=None)