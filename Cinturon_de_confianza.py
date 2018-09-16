# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 17:47:45 2018

@author: Samantha
"""

#codigo en elaboración, ya va


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



#Cantidad de valores que voy a generar en cada Montecarlo
cant = int(1e4)

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

# scipy.stats.uniform devuelve uniforme entre loc y loc+scale
Rm = c14m_sim/uniform.rvs(loc=c12im_p, scale=c12fm_p-c12im_p, size=cant)
Rf = c14f_sim/uniform.rvs(loc=c12if_p, scale=c12ff_p-c12if_p, size=cant)
# Para el estándar, la corriente final es más chica que la inicial
# Por lo tanto los extremos de la uniforme son al revés que en los otros casos
Rstd = c14std_sim/uniform.rvs(loc=c12fstd_p, scale=c12istd_p-c12fstd_p, size=cant)

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



edad = -tau*np.log(Rtot)
# Construyo el histograma de la variable Edad
Emin, Emax = np.min(edad), np.max(edad)

puntos_b = np.linspace(Emin, Emax,300) #puntos para graficar f(x)
bines = np.linspace(Emin, Emax,100)
numero, bins = np.histogram(edad, bins = bines) #numero=numero de entradas por bin
error = np.sqrt(numero) / (np.diff(bins)* np.sum(numero)) #error poissoniano
numero = numero / (np.diff(bins) * np.sum(numero)) #Normalizo a 1 (divido por el área ocupada por el histograma)

fig = plt.figure(figsize=(10,6))
plt.bar(bins[:-1], numero, width = np.diff(bins), yerr = error, ecolor="b", color='c', alpha=0.7)
#plt.plot(puntos_b, beta.pdf(puntos_b, ajuste_t[0], ajuste_t[1], loc=ajuste_t[2], scale=ajuste_t[3]), 'm-')
#plt.legend(loc=1, borderaxespad=0.)
plt.xlim([Emin, Emax])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Edad (años)', fontsize=12)
plt.ylabel('P(Edad)', fontsize=12)
#plt.title('Histograma')
plt.grid()
plt.show()



def esp(b,n):
    return np.sum(np.diff(b)*n*b[:-1])

esperanza = esp(bines, numero)

def lado_izq(hasta, bines, numero, bins):
    cerca = np.argmax([x for x in bines if x < hasta])
    return np.sum(np.diff(bins[:cerca]) * (numero[:cerca-1]))

def lado_der(desde, bines, numero, bins):
    cerca = np.argmax([x for x in bines if x >= desde])
    return np.sum(np.diff(bins[cerca:]) * (numero[cerca:]))

#%%
R_medido = esperanza

edades, probabilidad = [], []
for i in range(0,200):
    c14m_sim = poisson.rvs(np.sum(c14m)-1000+10*i,size=cant)/np.array(len(c14m), float)
    Rm = c14m_sim/uniform.rvs(loc=c12im_p, scale=c12fm_p-c12im_p, size=cant)
    Rtot = (Rm-Rf)/(Rstd-Rf)
    mask = (Rtot != np.float('+inf'))
    Rtot = Rtot[mask]
    mask = (Rtot>=0.00)
    Rtot = Rtot[mask]
    edad = -tau*np.log(Rtot)
    Emin, Emax = np.min(edad), np.max(edad)
    bines = np.linspace(Emin, Emax,100)
    numero, bins = np.histogram(edad, bins = bines)
    numero = numero / (np.diff(bins) * np.sum(numero))
    esperanza = esp(bines, numero)
    lado01 = lado_der(R_medido, bines, numero, bins)
    edades.append(int(np.mean(edad)))
    probabilidad.append(lado01)
    print(i, int(np.mean(edad)), int(esperanza), round(lado01,2))
#%%
    
#Esto a veces flashea, conviene mirar la lista
def find_index_of_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

indice_16 = find_index_of_nearest(probabilidad, 0.16)
indice_84 = find_index_of_nearest(probabilidad, 0.84)
print('edad_16', edades[indice_16])
print('edad_84', edades[indice_84])
