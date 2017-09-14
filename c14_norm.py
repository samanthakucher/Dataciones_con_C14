import numpy as np
from scipy.stats import poisson, norm
import matplotlib.pyplot as plt

num_muestras = int(1e4)
mu_muestra, mu_std = 212, 178
c14_muestra = poisson.rvs(mu_muestra, size=num_muestras)
c14_std = poisson.rvs(mu_std, size=num_muestras)
c12_muestra, c12_std = 2.75e12, 1e12
R_muestra = c14_muestra/c12_muestra
R_std = c14_std/c12_std
    
R = np.array(R_muestra, dtype=float)/np.array(R_std,dtype=float)

#Elimino los infinitos
mask = (R != np.float('+inf'))
R = R[mask]

# Verifico que para todos los elementos sea cierto que no son +inf
np.all(R != np.float('+inf'))

#Grafico la distribucion y la aproximo por una gaussiana
mean,std=norm.fit(R)
'''
fig = plt.figure(figsize=(20,10))
plt.hist(R, bins=60, normed=True, color='c')
plt.title('Relacion isotopica')
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
y = norm.pdf(x, mean, std)
plt.plot(x, y, 'm-', lw=2)
plt.grid()
plt.show()
'''
intervalo = norm.interval(alpha=0.68, loc=mean, scale=std)

tau = 5730./np.log(2)

def edad(rt):
    return -tau*np.log(rt)

edades = edad(R)
fig = plt.figure(figsize=(20,10))
mean,std=norm.fit(edades)
intervalo_edades = norm.interval(alpha=0.68, loc=mean, scale=std)
plt.hist(edades, bins=60, normed=True, color='c')
plt.title('Edades')
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
y = norm.pdf(x, mean, std)
plt.plot(x, y, 'm-', lw=2)
plt.plot(x, norm.pdf(x, intervalo_edades[0], std), 'g-')
plt.plot(x, norm.pdf(x, intervalo_edades[1], std), 'g-')
plt.grid()
plt.show()
intervalo_edades = norm.interval(alpha=0.68, loc=mean, scale=std)

R_medida = (mu_muestra/c12_muestra)/(mu_std/c12_std)
t_medida = edad(R_medida)
# Como la funcion tiene un -, el R_max se traduce en el t_min y viceversa
t_max  = edad(intervalo[0])
t_min = edad(intervalo[1])
print('Edad medida = {}, edad minima = {}, edad maxima = {}').format(t_medida, intervalo_edades[0], intervalo_edades[1])
print('Error superior = {}, error inferior = {}').format(t_medida-intervalo_edades[0], intervalo_edades[1]-t_medida)