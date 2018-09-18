# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:53:45 2018

@author: Samantha
"""
# Este código calcula los conteos de 12C y 14C obtenidos para la muestra
# con el estandar y con el fondo y a partir de ellos la relación isotópica
# en la muestra. Para hacerlo usa datos de mediciones reales cargados al inicio.

# Luego hace montecarlo para hallar la distribución de la variable aleatoria Rtot que
# resulta de dividir el cociente de la poissoniana del conteo de 14C en la muestra
# por la uniforme de la corriente de 12C con el mismo cociente para el estandar.

# Por último ajusta una pdf beta sobre la distribución y grafica el histograma 
# con el ajuste superpuesto.

import numpy as np
from scipy.stats import poisson, uniform, beta, gamma
import matplotlib.pyplot as plt
import time

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

# Lo que sigue es válido sólo cuando las mediciones son todas de igual duración.
# Aunque no es necesario tener el mismo número de mediciones para todos los casos.

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

#Cantidad de valores que voy a generar en cada Montecarlo
cant = int(1e7)

start_time = time.time()

# Asumimos conteo poissoniano 
# y que la corriente tiene distribución uniforme entre el valor inicial y final
# se toma el promedio como mejor estimador del u de la poissoniana.
# Asi como esta planteado el promedio también será la varianza, pero si 
# x tiene distribución de Poisson, que distribución tendrá su promedio?

#La suma de variables aleatorias poissonianas es poissoniana con parametro \sum \lambda_i
c14m_sim = poisson.rvs(np.sum(c14m),size=cant)/np.array(len(c14m), float)
c14f_sim = poisson.rvs(np.sum(c14f),size=cant)/np.array(len(c14f), float)
c14std_sim = poisson.rvs(np.sum(c14std),size=cant)/np.array(len(c14std), float)

# Acá ponemos los C12 fijos
Rm = c14m_sim / ((c12im_p + c12fm_p) / 2)
Rf = c14f_sim/ ((c12if_p + c12ff_p) / 2)
# Para el estándar, la corriente final es más chica que la inicial
# Por lo tanto los extremos de la uniforme son al revés que en los otros casos
Rstd = c14std_sim / ((c12fstd_p + c12istd_p) / 2)

Rtot = (Rm-Rf)/(Rstd-Rf)
# Lo correcto sería restar conteos con fondos primero y luego dividir.
# Para hacerlo hay que tener tener conteos y fondos por unidad de tiempo y de corriente. Por?
# De todos modos sugiero dejarlo así y discutirlo después. Ok

#Elimino los infinitos
mask = (Rtot != np.float('+inf'))
Rtot = Rtot[mask]
# Elimino los posibles casos en que esta variable tome valores negativos
mask = (Rtot>=0.00)
Rtot = Rtot[mask]

# Verifico que para todos los elementos sea cierto que no son +inf
np.all(Rtot != np.float('+inf'))
# Verifico que para todos los elementos sea cierto que no son negativos
np.all(Rtot <= 0.00)

print("--- %s seconds ---" % (time.time() - start_time))
#%%

# Construyo el histograma de la variable Rtot
Rmin, Rmax = 0, 1.25

puntos_a = np.linspace(Rmin, Rmax,300) #puntos para graficar f(x)
bines = np.linspace(Rmin, Rmax,100)
numero, bins = np.histogram(Rtot, bins = bines) #numero=numero de entradas por bin
error = np.sqrt(numero) / (np.diff(bins)* np.sum(numero)) #error poissoniano
numero = numero / (np.diff(bins) * np.sum(numero)) #Normalizo a 1 (divido por el área ocupada por el histograma)

# Ajusto el histograma con una pdf beta
#ajuste = beta.fit(Rtot)

#%%
# Grafico el histograma con el ajuste superpuesto
fig = plt.figure(figsize=(10,6))
plt.bar(bins[:-1], numero, width = np.diff(bins), yerr = error, ecolor="b", color='c', alpha=0.7)
#plt.plot(puntos_a, beta.pdf(puntos_a, ajuste[0], ajuste[1], loc=ajuste[2], scale=ajuste[3]), 'm-')
#plt.legend(loc=1, borderaxespad=0.)
plt.xlim([Rmin, Rmax])
plt.xlabel('R')
plt.title('Histograma')
plt.grid()
plt.show()

#%%
# Ahora que se cuenta con una distribución para la variable R_tot hay que encontrar
# la distribución de la variable edad de la muestra asociada a esta variable R_tot.
# Como se conoce la distribución de R_tot y su relación analítica con la edad de la muestra
# la distribución de esta última puede escribirse analíticamente usando un cambio de variables.

# También se puede continuar con el Montecarlo y ver que distribución toma la variable 
# edad de la muestra ingresando R_tot random con la pdf beta hallada en la fórmula para
# calcular la edad de la muestra dado R_tot.

# Luego, se puede calcular la dispersión de la variable edad de la muestra y se tendrá su error.

#%%

# 1) Tomo el intervalo de 68% CL para Rtot y lo traduzco en un intervalo para la edad
# basandome en que la distribución de Rtot se ajusta bien por una beta

intervalo_r = beta.interval(0.68, ajuste[0], ajuste[1], loc=ajuste[2], scale=ajuste[3])
intervalo_edad_r = [-tau*np.log(intervalo_r[1]), -tau*np.log(intervalo_r[0])]

#%%

# 2) Veo qué distribucion toma la variable edad 

edad = -tau*np.log(Rtot)
# Construyo el histograma de la variable Edad
Emin, Emax = np.min(edad), np.max(edad)

puntos_b = np.linspace(Emin, Emax,300) #puntos para graficar f(x)
bines = np.linspace(Emin, Emax, 100)
numero, bins = np.histogram(edad, bins = bines) #numero=numero de entradas por bin
error = np.sqrt(numero) / (np.diff(bins)* np.sum(numero)) #error poissoniano
numero = numero / (np.diff(bins) * np.sum(numero)) #Normalizo a 1 (divido por el área ocupada por el histograma)

# Ajusto el histograma con una pdf beta
#ajuste_t = beta.fit(edad)



#%%

# Tomo el intervalo basandome en que edad tiene distribucion beta

intervalo_edad_hist = beta.interval(0.68, ajuste_t[0], ajuste_t[1], loc=ajuste_t[2], scale=ajuste_t[3])
#%%
# Tomo el intervalo numéricamente

intervalo_edad_num = np.array([np.percentile(edad, 16), np.percentile(edad, 84)], dtype=int)
# Grafico el histograma con las barras de los percentiles
bines = np.linspace(Emin, Emax, 1000)
numero, bins = np.histogram(edad, bins = bines) #numero=numero de entradas por bin
numero = numero / (np.diff(bins) * np.sum(numero)) #Normalizo a 1 (divido por el área ocupada por el histograma)

fig = plt.figure(figsize=(10,6))
plt.bar(bins[:-1], numero, width = np.diff(bins), ecolor="b", color='c', alpha=0.7)
plt.xlim([Emin, Emax])
plt.tick_params(labelsize=16)
plt.xlabel('Edad (años)', fontsize=18)
plt.ylabel('Densidad de eventos (adim.)', fontsize=18)
altura_lineas = [numero[np.where(bins >= intervalo_edad_num[i])[0][0] - 1] for i in [0,1]]
plt.vlines(intervalo_edad_num[0], ymin=0, ymax=altura_lineas[0], color='r')
plt.vlines(intervalo_edad_num[1], ymin=0, ymax=altura_lineas[1], color='r')


from matplotlib.ticker import FormatStrFormatter
ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0e'))

#plt.title('Histograma')
plt.grid()
plt.show()


#%%
# Grafico el histograma con el ajuste superpuesto
fig = plt.figure(figsize=(10,6))
plt.bar(bins[:-1], numero, width = np.diff(bins), yerr = error, ecolor="b", color='c', alpha=0.7)
#plt.plot(puntos_b, beta.pdf(puntos_b, ajuste_t[0], ajuste_t[1], loc=ajuste_t[2], scale=ajuste_t[3]), 'm-')
#plt.legend(loc=1, borderaxespad=0.)
plt.xlim([Emin, Emax])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Edad (años)', fontsize=12)
plt.ylabel('P(Edad)', fontsize=12)
plt.vlines(intervalo_edad_hist[0], ymin=0, ymax=beta.pdf(intervalo_edad_hist[0], ajuste_t[0], ajuste_t[1], loc=ajuste_t[2], scale=ajuste_t[3]), color='r')
plt.vlines(intervalo_edad_hist[1], ymin=0, ymax=beta.pdf(intervalo_edad_hist[1], ajuste_t[0], ajuste_t[1], loc=ajuste_t[2], scale=ajuste_t[3]), color='r')
#plt.title('Histograma')
plt.grid()
plt.show()

#%%

# Comparo ambos metodos
print('Intervalos para la edad calculados a partir de la distribucion de')

print('R = ', np.array(intervalo_edad_r,dtype=int))
vmedio = (intervalo_edad_r[0] + intervalo_edad_r[1]) / 2
error = np.abs(intervalo_edad_r[0] - intervalo_edad_r[1]) / 2
print('Error relativo {}%.'.format(error / vmedio))

print('la edad = ', np.array(intervalo_edad_hist,dtype=int))
vmedio = (intervalo_edad_hist[0] + intervalo_edad_hist[1]) / 2
error = np.abs(intervalo_edad_hist[0] - intervalo_edad_hist[1]) / 2
print('Error relativo {}%.'.format(error / vmedio))
# Ahora ya no da enorme! El intervalo es aproximadamente [8150 años, 8810 años]

print('percentiles = ', intervalo_edad_num)
vmedio = (intervalo_edad_num[0] + intervalo_edad_num[1]) / 2
error = np.abs(intervalo_edad_num[0] - intervalo_edad_num[1]) / 2
print('Error relativo {}%.'.format(error / vmedio))