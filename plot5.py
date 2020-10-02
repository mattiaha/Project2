# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
#declare the size of the matrix we have used
N = 4
#write in the name of the file here, from same folder as code is saved in

data1 = np.loadtxt("100N_QM2_r10.txt", skiprows=2)
data2 = np.loadtxt("100N_QM2_r50.txt", skiprows=2)
data3 = np.loadtxt("100N_QM2_r60.txt", skiprows=2)
data4 = np.loadtxt("100N_QM2_r70.txt", skiprows=2)

x = data1[:,0]
plt.plot(x,data1[:,4])

plt.plot(x,data2[:,4])
plt.plot(x,data3[:,4])
plt.plot(x,data4[:,4])
plt.xlabel("$\omega_{r}$")
plt.ylabel("Relative error")
#write in size of N
plt.title("Relative error for groundstates")
plt.legend(["$\\rho$ =10","$\\rho=50$","$\\rho$=60","$\\rho$=70"])

plt.show()