# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 16:50:52 2020

@author: matti
"""

import os
import matplotlib.pyplot as plt
import numpy as np
#declare the size of the matrix we have used
N = 100
#write in the name of the file here, from same folder as code is saved in

data1 = np.loadtxt("100N_QM_1.txt", skiprows=2)
data2 = np.loadtxt("100N_QM_5.txt", skiprows=2)
data3 = np.loadtxt("100N_QM_10.txt", skiprows=2)
data4 = np.loadtxt("100N_QM_15.txt", skiprows=2)

x = np.linspace(1,N,N)
plt.plot(x,data1[:,3])

plt.plot(x,data2[:,3])
plt.plot(x,data3[:,3])
plt.plot(x,data4[:,3])
plt.xlabel("Number of eigenvalue")
plt.ylabel("Relative error")
#write in size of N
plt.title("N=100")
plt.legend(["rho =1","rho=5","rho=10","rho=15"])

plt.show()