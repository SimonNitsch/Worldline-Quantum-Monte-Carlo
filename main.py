import numpy as np
import matplotlib.pyplot as plt
import sys
import os

arguments = sys.argv
foldername = arguments[1]

folder_path = os.path.join(os.getcwd(), foldername)
os.chdir(folder_path)
E_data = np.fromfile("E_data", np.float64)
E_data = E_data.reshape(-1, 2)
S_data = np.fromfile("S_data", np.float64)
S_data = S_data.reshape(-1, 2)
E = E_data[:,0]
Delta_E = E_data[:,1]
S = S_data[:,0]
Delta_S = S_data[:,1]


plt.figure()
plt.subplot(2,1,1)
plt.errorbar(np.arange(E.size),E,Delta_E)
plt.title("Energy over time")
plt.subplot(2,1,2)
plt.errorbar(np.arange(S.size),S,Delta_S)
plt.title("Spin-Spin-Correlation over time")
plt.show()

