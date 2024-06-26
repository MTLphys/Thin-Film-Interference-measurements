import matplotlib.pyplot as plt 
import numpy as np 

nkAlq   =  np.loadtxt('nk_Alq3.txt',skiprows =0) 
nptcda1 =  np.loadtxt('n_PTCDA_Niranjala (2 par and 2 per).txt',skiprows =1)
nptcda2 =  np.loadtxt('PTCDA (FDTD) nkparSi.txt',skiprows =0)
nptcda3 =  np.loadtxt('PTCDA (FDTD).txt',skiprows =0)



plt.plot(nkAlq[:,0],nkAlq[:,1],label = 'Alq3 n')
plt.plot(nkAlq[:,0],nkAlq[:,2],label = 'Alq3 k')
plt.plot(nptcda1[:,0],nptcda1[:,1],label = 'PTCDA n1a')
plt.plot(nptcda1[:,0],nptcda1[:,2],label = 'PTCDA n2a')
plt.plot(nptcda1[:,0],nptcda1[:,3],label = 'PTCDA k1a')
plt.plot(nptcda1[:,0],nptcda1[:,4],label = 'PTCDA k2a')

plt.plot(nptcda2[:,0],nptcda2[:,1],label = 'PTCDA nb')
plt.plot(nptcda2[:,0],nptcda2[:,2],label = 'PTCDA kb')

plt.plot(nptcda3[:,0],nptcda3[:,1],label = 'PTCDA nc')
plt.plot(nptcda3[:,0],nptcda3[:,2],label = 'PTCDA kc')

plt.legend()

plt.show()