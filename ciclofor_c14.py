# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

# A veces este programa tira error, hay que ejecutarlo hasta que ande

num_muestras = int(1e3) #cantidad de variables aleatorias que voy a generar
mu_std = 178 
num_mediciones = 10
c12_muestra, c12_std = 2.75e12, 1e12 # No son variables aleatorias
R_std_prom = mu_std/c12_std
valor_actual = np.mean(R_std_prom)*c12_muestra
edadmax_gauss, edadmin_gauss = [], []
edadmax_sim, edadmin_sim = [], []
edad_gauss, edad_sim = [], []
for mu_muestra in range(1, int(valor_actual)):
    c14_muestra = poisson.rvs(mu_muestra, size=num_muestras)
    mean14,std14 = norm.fit(c14_muestra) # Aproximo por una gaussiana
    # Sup que los mu de las muestras tienen esta distribución

    c14_std = poisson.rvs(mu_std, size=num_muestras)
    mean14std,std14std = norm.fit(c14_std) # Aproximo por una gaussiana
    # Sup que los mu de las std tienen esta distribución
    R_muestra = c14_muestra/c12_muestra
    R_std = c14_std/c12_std
    R = np.array(R_muestra, dtype=float)/np.array(R_std,dtype=float)
    #Elimino los infinitos
    mask = (R != np.float('+inf'))
    R = R[mask]
    # Verifico que para todos los elementos sea cierto que no son +inf
    np.all(R != np.float('+inf'))

    gaussR, sigmagaussR = norm.fit(R) #Esto es trabajando con un unico valor
    #Ahora simulamos todas las mediciones
    mediciones = np.zeros((num_mediciones, num_muestras))
    for i in range(0, num_mediciones):
        # Esto simula cada medicion por separado
        mu_medicion_muestra = int(norm.rvs(loc=mean14, scale=std14)) 
        mu_medicion_std = int(norm.rvs(loc=mean14std, scale=std14std))
        medicion_muestra = poisson.rvs(mu_medicion_muestra, size=num_muestras)
        medicion_std = poisson.rvs(mu_medicion_std, size=num_muestras)
        R_muestra = medicion_muestra/c12_muestra
        R_std = medicion_std/c12_std
        B = np.array(R_muestra, dtype=float)/np.array(R_std,dtype=float)
        mask = (B != np.float('+inf')) #Elimino los infinitos
        B = B[mask]
        np.all(B != np.float('+inf')) # Verifico que sea cierto que no son +inf
        mediciones[i,:] = B
    # Cada fila de la matriz mediciones es una medición. Las promedio:
    def prom_matriz(A, f, c):
        suma = np.zeros(c)
        for i in range(0,f):
            suma = A[i,:]+suma
        return suma/f
    
    Bprom = prom_matriz(mediciones, num_mediciones, num_muestras)
    gaussprom, sigmaprom = norm.fit(Bprom) #Aproximo este promedio por una gaussiana
    intervalo_gauss = norm.interval(alpha=0.68, loc=gaussR, scale=sigmagaussR/np.sqrt(num_mediciones))
    intervalo_sim = norm.interval(alpha=0.68, loc=gaussprom, scale=sigmaprom)

    #Traduzco estos valores a edad:
    tau = 5730./np.log(2)
    def edad(rt):
        return -tau*np.log(rt)
    
    # Como la funcion tiene un -, el R_max se traduce en el t_min y viceversa
    t_gauss, t_sim = edad(gaussR), edad(gaussprom)
    t_max_gauss, t_min_gauss  = edad(intervalo_gauss[0]), edad(intervalo_gauss[1])
    t_max_sim, t_min_sim  = edad(intervalo_sim[0]), edad(intervalo_sim[1])
    edad_gauss.append(t_gauss)
    edad_sim.append(t_sim)
    edadmax_gauss.append(t_max_gauss)
    edadmin_gauss.append(t_min_gauss)
    edadmax_sim.append(t_max_sim)
    edadmin_sim.append(t_min_sim)

#%%

fig = plt.figure()
plt.plot(edad_gauss, edadmax_gauss, 'm-')
plt.plot(edad_gauss, edadmin_gauss, 'm-', label='Intervalo')
plt.plot(edad_gauss, edad_gauss, 'k-', label='Edad')
#plt.plot(edad_sim, edad_sim, 'k-')
#plt.plot(edad_sim, edadmax_sim, 'g-')
#plt.plot(edad_sim, edadmin_sim, 'g-', label='Intervalo')
plt.grid()
plt.xlabel('Edad')
plt.ylabel('Intervalo 68% CL')
plt.legend(loc=0)
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.show()

#plt.savefig('F:/Facultad/Estadística/ciclofor.png', dpi=None, facecolor='w', edgecolor='w',         orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1, frameon=None)

#%%
delta_gauss = (np.array(edadmax_gauss)-np.array(edadmin_gauss))/2.
delta_sim = (np.array(edadmax_sim)-np.array(edadmin_sim))/2.

fig = plt.figure()
plt.plot(edad_gauss, delta_gauss, 'm-')
#plt.plot(edad_sim, delta_sim, 'g-')
plt.xlim(xmin=0)
plt.grid()
plt.xlabel('Edad')
plt.ylabel('$\Delta$Edad')
plt.show()

#plt.savefig('F:/Facultad/Estadística/ciclofor_delta.png', dpi=None, facecolor='w', edgecolor='w',         orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1, frameon=None)