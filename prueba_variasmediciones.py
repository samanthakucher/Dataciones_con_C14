import numpy as np
from scipy.stats import poisson, norm
import random as rd
import matplotlib.pyplot as plt

# La idea es simular varias mediciones (si las tienen las pueden ecribir directamente)
# Como no las tengo, supongo que los mus de las poisson vienen de unas gaussianas 
# centradas en el mu_muestra y el mu_std 

num_muestras = int(1e4) #para simular
mu_muestra, mu_std = 2, 178 #promedio 212 muestra
num_mediciones = 20
c14_muestra = poisson.rvs(mu_muestra, size=num_muestras)

mean14,std14 = norm.fit(c14_muestra) 
intervalo_c14 = norm.interval(alpha=0.68, loc = mean14, scale=std14)
# Sup que los lambda de la muestra se mueven en este intervalo

c14_std = poisson.rvs(mu_std, size=num_muestras)
mean14std,std14std = norm.fit(c14_std) 
intervalo_c14std = norm.interval(alpha=0.68, loc = mean14std, scale=std14std)
# Sup que los lambda del std se mueven en este intervalo

c12_muestra, c12_std = 2.75e12, 1e12 # No son variables aleatorias


R_muestra = c14_muestra/c12_muestra
R_std = c14_std/c12_std
R = np.array(R_muestra, dtype=float)/np.array(R_std,dtype=float)
#Elimino los infinitos
mask = (R != np.float('+inf'))
R = R[mask]
# Verifico que para todos los elementos sea cierto que no son +inf
np.all(R != np.float('+inf'))

gaussR, sigmagaussR = norm.fit(R)

mediciones = np.zeros((num_mediciones, num_muestras))
for i in range(0, num_mediciones):
    # Esto simula cada medicion por separado
    mu_medicion_muestra = rd.randint(int(intervalo_c14[0]), int(intervalo_c14[1]))
    mu_medicion_std = rd.randint(int(intervalo_c14std[0]), int(intervalo_c14std[1]))
    medicion_muestra = poisson.rvs(mu_medicion_muestra, size=num_muestras)
    medicion_std = poisson.rvs(mu_medicion_std, size=num_muestras)
    R_muestra = medicion_muestra/c12_muestra
    R_std = medicion_std/c12_std
    B = np.array(R_muestra, dtype=float)/np.array(R_std,dtype=float)
    mask = (B != np.float('+inf')) #Elimino los infinitos
    B = B[mask]
    np.all(B != np.float('+inf')) # Verifico que sea cierto que no son +inf
    mediciones[i,:] = B
    
#gaussB, sigmagaussB = norm.fit(B)    

#Bprom = np.sum(mediciones[i,:])/num_mediciones
def prom_matriz(A, f, c):
    suma = np.zeros(c)
    for i in range(0,f):
        suma = A[i,:]+suma
    return suma/f
    
Bprom = prom_matriz(mediciones, num_mediciones, num_muestras)
gaussprom, sigmaprom = norm.fit(Bprom)

x = np.linspace(np.amin(R),np.amax(R),200)
fig = plt.figure(figsize=(20,10))
plt.hist(R, bins=60, normed=True, color='m', alpha=0.5)
#plt.hist(B, bins=60, normed=True, color='c', alpha = 0.5)
plt.hist(Bprom, bins=60, normed=True, color='y', alpha = 0.5)
plt.plot(x, norm.pdf(x, gaussprom, sigmaprom), 'g-')
plt.plot(x, norm.pdf(x, gaussR, sigmagaussR/np.sqrt(num_mediciones)), 'b-')
plt.plot(x, norm.pdf(x, gaussR, sigmagaussR), 'r-')
#plt.plot(x, norm.pdf(x, gauss2, sigmagauss2), 'b-')
plt.grid()
plt.show()

print(gaussR, sigmagaussR)
#print(gauss2, sigmagauss2)
print(gaussprom, sigmaprom)
