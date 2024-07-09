

import matplotlib.pyplot as plt 
import numpy as np 
"""
section for loading and setting parameters of the fitting 

"""
#applying a silly thing


nkAlq   =  np.loadtxt('nk_Alq3.txt',skiprows =0) 
nptcda1 =  np.loadtxt('n_PTCDA_Niranjala (2 par and 2 per).txt',skiprows =1)
nptcda2 =  np.loadtxt('PTCDA (FDTD) nkparSi.txt',skiprows =0)
nptcda3 =  np.loadtxt('PTCDA (FDTD).txt',skiprows =0)
nptcdaavg =(nptcda1[:,1]+nptcda1[:,2])/2
nptcda4=np.loadtxt('PTCDAInterpolatedn.txt',skiprows =0)
nbk7 =  np.loadtxt('N-BK7.txt',skiprows =1)

d1g = 300 #fit parameters 1 for layer thickness 
d2g = 000
d3g = 0 


"""
section for interpolation of parameters

"""
def nkinterp(ndata,kdata,wldata,wl,wlkdata=np.zeros(1)):
    """_linearly interpolate data for the complex and real index 
    of refraction

    Args:
        ndata (array):  the index of refraction from data
        kdata (array):  the extinction coef from data 
        wldata (array): the wavelength from the data
        wl (array): the target wavelength space
        wlkdata (array, optional): defaults to wldata if undefined the data for the k 

    Returns:
        _type_: interpolated index of refraction
    """
    if(len(wlkdata)==1):
        print('all 0')
        wlkdata = wldata
    ntil = np.interp(wl,wldata,ndata)+1.0j*np.interp(wl,wlkdata,kdata)
    return ntil
    
wl = np.linspace(200,900,3000)
print(nptcda1.shape)
n0 = 1.0004*wl/wl
n1 = nkinterp(nptcda4[:,1],nptcda4[:,2],nptcda4[:,0],wl)
n2 = nkinterp(nkAlq[:,1],nkAlq[:,2],nkAlq[:,0],wl)
#n3=  nkinterp(nptcda2[:,1],nptcda2[:,2],nptcda2[:,0],wl)
n3 = nkinterp(nbk7[0:101,1],nbk7[102:125,1],nbk7[0:101,0]*1000,wl,nbk7[102:125,0]*1000)#1.45*np.ones_like(n2)#
"""
Section for showing the plots of the indeces of refraction


"""

shownkplots =  False# True

if(shownkplots):
    fig,ax = plt.subplots()
    ax.plot(nkAlq[:,0],nkAlq[:,1],label = 'Alq3 n')
    ax.plot(nkAlq[:,0],nkAlq[:,2],label = 'Alq3 k')
    ax.plot(wl,np.real(n1),label='Alq3 n interp')
    ax.plot(wl,np.imag(n1),label='Alq3 k interp')
    #ax.plot(nptcda1[:,0],nptcda1[:,1],label = 'PTCDA n1a')
    #ax.plot(nptcda1[:,0],nptcda1[:,2],label = 'PTCDA n2a')
    #ax.plot(nptcda1[:,0],nptcda1[:,3],label = 'PTCDA k1a')
    #ax.plot(nptcda1[:,0],nptcda1[:,4],label = 'PTCDA k2a')

    ax.plot(nptcda2[:,0],nptcda2[:,1],label = 'PTCDA nb')
    ax.plot(nptcda2[:,0],nptcda2[:,2],label = 'PTCDA kb')
    ax.plot(wl,np.real(n2),label='PTCDA n interp')
    ax.plot(wl,np.imag(n2),label='PTCDA k interp')
    #ax.plot(nptcda3[:,0],nptcda3[:,1],label = 'PTCDA nc')
    #ax.plot(nptcda3[:,0],nptcda3[:,2],label = 'PTCDA kc')
    ax.plot(1000*nbk7[0:101,0],nbk7[0:101,1],label = 'Glass nb')
    ax.plot(1000*nbk7[102:125,0],nbk7[102:125,1],label = 'Glass kb')
    ax.plot(wl,np.real(n3),label='Glass n interp')
    ax.plot(wl,np.imag(n3),label='Glass k interp')

    ax.legend()


"""section for calculating the M Matrices given above parameters"""
def rij(ni,nj):
    """The Fresnell reflection coeffiecent

    Args:
        ni (complex_128): index of refraction for layer i
        nj (complex_128): index of refraction for layer j 

    Returns:
        complex_128: reflection coeffiecent matrix
    """
    return (ni-nj)/(ni+nj)

def tij(ni,nj):
    """The fresnel transmission coeffeicient 

    Args:
        ni (complex_128):index of refraction for layer i
        nj (complex_128):index of refraction for layer i

    Returns:
        _type_: transmission coeffiecent matrix
    """
    return 2*(ni)/(ni+nj)
def delta(ni,di,wl):
    """single layer phase delay

    Args:
        ni (_type_): index of refraction 
        di (_type_): thickness of layer 
        wl (_type_): wavelength
    Returns: 
        complex_128: the phase difference for a single layer i 
    """
    return np.pi*2*ni*di/wl
def dot(A,B):
    C=np.zeros_like(A)
    for i in range(A.shape[2]):
        C[:,:,i]=A[:,:,i].dot(B[:,:,i])
    return C
def Mij(ni,nj,dj,wl):
    """Abbles Matrix for a single layer of material 

    Args:
        tij (complex 128): _description_
        rij (_type_): _description_
        deltaj (_type_): _description_

    Returns:
        _type_: _description_
    """
    return 1/tij(ni,nj)*np.asarray([[           np.exp(-1j*delta(nj,dj,wl)),rij(ni,nj)*np.exp(1j*delta(nj,dj,wl)) ],
                                    [rij(ni,nj)*np.exp(-1j*delta(nj,dj,wl)),           np.exp(1j*delta(nj,dj,wl)) ]])


def M(d,n,wl,layers):
    Mtot = np.zeros((2,2,len(wl)))
    for i in range(len(layers)-2):
        
        if i==0 :
            Mtot = Mij(n[layers[i]],n[layers[i+1]],d[i+1],wl)
        else:
            Mtot = dot(Mtot,Mij(n[layers[i]],n[layers[i+1]],d[i+1],wl))
    rback = rij(n[layers[len(layers)-2]],n[layers[len(layers)-1]])
    rbmat = np.asarray([[np.ones_like(rback),rback],
                      [rback,np.ones_like(rback)]])
    Mtot = dot(Mtot,rbmat)
    return 1/tij(n[layers[len(layers)-2]],n[layers[len(layers)-1]])*Mtot

layers = [0,2, 2, 3]     #the material type of each sequential layer
d =      [0,d1g,d1g,  0] #the thickness of each layer, start and end layers set to 0 as no transfer matrix is applied
n =      [n0,n1, n2, n3] #air,ptcda,alq3,glass

#print(M(d,n,wl,layers))
#delta1 = delta(n[2],d[1],wl)
# delta2 = delta(n[2],d[2],wl)
# delta3 = delta(n[1],d[3],wl)
# M1 = Mij(n[0],n[2],d[1],wl)
# # M2 = Mij(n[1],n[2],d[2],wl)
# # M3 = Mij(n[2],n[1],d[3],wl)
# t34 = tij(n[2],n[3])
# r34 = rij(n[2],n[3])
# mr34 = np.asarray([[np.ones_like(r34),r34],
#                    [r34,np.ones_like(r34)]])
# Mfa = 1/t34*dot(dot(dot(M1,M2),M3),mr34)
Mfb = M(d,n,wl,layers)
fig, ax = plt.subplots()
#ax.plot(wl,np.abs(Mfa[1,0,:]/Mfa[0,0,:])**2)
#ax.plot(wl,np.abs(1/Mfa[0,0,:])**2)

ax.plot(wl,np.abs(Mfb[1,0,:]/Mfb[0,0,:])**2/.041,label ="Reflectance" )
ax.plot(wl,np.abs(1/Mfb[0,0,:])**2/.64, label ="Transmittance" )


#ax.plot(wl,np.abs(r34))
#ax.plot(wl,np.abs(t34))
ax.set_ylim(-.1,1.2)
#ax.set_ylim(0,17)
ax.legend()
plt.show()