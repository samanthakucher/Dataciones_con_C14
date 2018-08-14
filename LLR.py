# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 21:31:06 2018

@author: Samantha
"""

import numpy as np
from scipy.stats import poisson, uniform, beta
from scipy.optimize import curve_fit
from scipy.special import gamma
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import time


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

# Hay 22 mediciones de c14m, 8 de c14std y 9 de c14f
# Supongo que "cada medicion" es de la muestra
# Entonces estimo Rstd y Rfondo directo con todas las mediciones juntas
# y a la likehood le voy agregando mediciones de la muestra

cant = int(1e4)
c12if_p, c12ff_p = np.mean(c12if), np.mean(c12ff)
c12istd_p, c12fstd_p = np.mean(c12istd), np.mean(c12fstd)

c14f_sim = poisson.rvs(np.sum(c14f),size=cant)/np.array(len(c14f), float)
c14std_sim = poisson.rvs(np.sum(c14std),size=cant)/np.array(len(c14std), float)
Rf = c14f_sim/uniform.rvs(loc=c12if_p, scale=c12ff_p-c12if_p, size=cant)
Rstd = c14std_sim/uniform.rvs(loc=c12fstd_p, scale=c12istd_p-c12fstd_p, size=cant)
#%%
# Voy simulando las mediciones de a una y ajusto por una beta, cuyos
# parametros guardo en 4 vectores
start_time = time.time()
R, param_beta_0, param_beta_1, param_beta_2, param_beta_3 = [],[],[],[],[]
for i in range(1,len(c14m)+1):
    Rm = poisson.rvs(np.mean(c14m[:i]), size=cant)/uniform.rvs(loc=min(np.mean(c12im[:i]), np.mean(c12fm[:i])), scale=abs(np.mean(c12fm[:i])-np.mean(c12im[:i])), size=cant)
    Rtot = (Rm-Rf)/(Rstd-Rf)
    #Elimino los infinitos
    mask = (Rtot != np.float('+inf'))
    Rtot = Rtot[mask]
    # Elimino los posibles casos en que esta variable tome valores negativos
    mask = (Rtot>=0.00)
    Rtot = Rtot[mask]
    R.append(np.mean(Rtot))
    ajuste = beta.fit(Rtot)
    param_beta_0.append(ajuste[0])
    param_beta_1.append(ajuste[1])
    param_beta_2.append(ajuste[2])
    param_beta_3.append(ajuste[3])
print("--- %s seconds ---" % (time.time() - start_time))
#%%

aa,bb = np.linspace(min(param_beta_0),max(param_beta_0),100), np.linspace(min(param_beta_1),max(param_beta_1), 100)
volumen = np.zeros(shape=(100,100, len(R)))
for k in range(len(R)):
    print('k=',k)
    for s in range(100):
        for j in range(100):
            volumen[i, j, k] = beta.pdf(R[k], aa[s], bb[j])
#%%
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
(x, y) = np.meshgrid(aa,bb)
for k in range(len(R)):
    ax.plot_surface(x,y,volumen[:,:,k], cmap=cm.brg)
    ax.set_title(i)
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    plt.grid()
    plt.draw()
    plt.pause(1)
    plt.cla()

#%%
k=0   
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x,y,volumen[:,:,k], cmap=cm.brg)