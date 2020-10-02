# -*- coding: utf-8 -*-


import os
import matplotlib.pyplot as plt
import numpy as np
#declare the size of the matrix we have used
N = 4
#write in the name of the file here, from same folder as code is saved in
data = np.loadtxt("100N_QM2_r50.txt", skiprows=2)
x  = data[:,0]
approx= data[:,1]
exact = data[:,2]


plt.plot(x,approx,'b')
plt.plot(x,exact,'r')
#plot of results from armadillo and jacobi method overlap
plt.plot(x,approx,'bo')
plt.plot(x,exact,'ro')
plt.xlabel("$\omega_{r}$")
plt.ylabel("Eigenvalue")
#write in size of N
plt.title("Groundstates, $\\rho$=50")
plt.legend(["jacobi","exact"])
plt.show()