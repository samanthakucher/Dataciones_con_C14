import numpy as np
'''
#a)
A = np.matrix('1 -0.6 0.36; 1 -0.2 0.04; 1 0.2 0.04; 1 0.6 0.36')
V = np.matrix('4 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 4')
At = A.getT()
Vi = V.getI()
mult = At*Vi*A
Va = mult.getI()
y = np.matrix('5 3 5 8')
yt = y.getT()
param = Va*At*Vi*yt
print(param)
'''
#b)
A = np.matrix('1 -0.6 0.36; 1 -0.2 0.04; 1 0.2 0.04; 1 0.6 0.36')
At = A.getT()
mult = At*A
Va = mult.getI()
y = np.matrix('5 3 5 8')
yt = y.getT()
param = Va*At*yt
s1 = np.array(yt[0] -A[0]*param)
s2 = np.array(yt[1]- A[1]*param)
s3 = np.array(yt[2]-A[2]*param)
sigma2 = s1**2 + s2**2 + s3**2
sigma = np.sqrt(sigma2)
print(sigma)



