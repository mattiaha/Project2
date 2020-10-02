# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:02:58 2020

@author: matti
"""

import os
import matplotlib.pyplot as plt
import numpy as np
#declare the size of the matrix we have used
N = 10
#write in the name of the file here, from same folder as code is saved in
data = np.loadtxt("10N_QM.txt", skiprows=2)
approx  = data[:,0]
exact= data[:,1]
arma = data[:,2]
rel_err = data[:,3]
x = np.linspace(1,N,N)
plt.plot(x,approx)
plt.plot(x,exact)
#plot of results from armadillo and jacobi method overlap
"""plt.plot(x,armadillo)"""
plt.xlabel("Number of eigenvalue")
plt.ylabel("Eigenvalue")
#write in size of N
plt.title("N=10")
plt.legend(["jacobi","exact"])
plt.show()