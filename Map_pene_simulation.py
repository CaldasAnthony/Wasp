import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import matplotlib.cm as cm
from pydataread import *
from pyfunction import *

#alpha = np.load('alpha_n.npy')
#data = open('I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.dat','r')
alpha = np.load('Files_penetration_duo/WASP121/alpha_H2O_CO2_var_2.npy')
data = open('WASP121_var_2.dat','r')

I = data.readlines()

wl_2 = np.zeros(len(I))
fl = np.zeros(len(I))
error = np.zeros(len(I))
for i in range(len(I)) :

    wl_2[i] = line_search(I[i])[0]
    fl[i] = line_search(I[i])[1]
wh_wl, = np.where((wl_2 > 1.)*(wl_2 < 10.))

bande_sample = np.load('/Users/caldas/Pytmosph3R/Source/bande_sample_bin10.npy')
wl = 10**4/bande_sample
obs_range, = np.where((wl > 1.)*(wl < 10.))

print obs_range.size
print wh_wl.size

wh_p = np.array([],dtype=np.int)
for i in range(wh_wl.size) :
    who, = np.where(np.round(wl_2[wh_wl[i]],5) == np.round(wl[obs_range],5))
    if who.size != 0 :
        wh_p = np.append(wh_p,np.array([who[0]]))

print wh_p.size

val = np.array([2, 3, 4, 5, 6, 7, 8, 9])
stick = np.array([])
for i in range(val.size) :
    wh, = np.where(wl[obs_range][::-1] >= val[i])
    stick = np.append(stick,np.array([wh[0]]))

print np.shape(alpha)
s = np.shape(alpha)
theta_res = s[2]

vin = np.amin(alpha)
vax = np.amax(alpha)

fact = theta_res/(np.amax(fl[wh_wl]) - np.amin(fl[wh_wl]))

plt.figure(1,figsize=(14,8))
plt.imshow(np.transpose(alpha[0,::-1,:]), aspect ='auto',vmin= vin,vmax=vax,cmap=plt.cm.get_cmap('nipy_spectral_r'))
plt.colorbar()
plt.plot(wh_p[::-1],(fl[wh_wl]-np.amin(fl[wh_wl]))*fact,'k',linewidth=2)
plt.yticks(np.linspace(0,theta_res,9))
plt.xticks(stick)
plt.ylim(0,theta_res-1)
plt.xlim(0,s[1])

plt.figure(2,figsize=(14,8))
plt.imshow(np.transpose(alpha[1,::-1,:]), aspect ='auto',vmin= vin,vmax=vax,cmap=plt.cm.get_cmap('nipy_spectral_r'))
plt.colorbar()
plt.plot(wh_p[::-1],(fl[wh_wl]-np.amin(fl[wh_wl]))*fact,'k',linewidth=2)
plt.yticks(np.linspace(0,theta_res,9))
plt.xticks(stick)
plt.ylim(0,theta_res-1)
plt.xlim(0,s[1])

vin = np.amin(np.transpose(alpha[0,::-1,:]) - np.transpose(alpha[1,::-1,:]))
vax = np.amax(np.transpose(alpha[0,::-1,:]) - np.transpose(alpha[1,::-1,:]))

'''
val = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9])
stick = np.array([])
for i in range(val.size) :
    wh, = np.where(wl[obs_range][::-1] >= val[i])
    stick = np.append(stick,np.array([wh[0]]))

plt.figure(3,figsize=(14,8))
plt.imshow(np.transpose(alpha[0,::-1,:][:500,:]) - np.transpose(alpha[1,::-1,:][:500,:]), aspect ='auto',vmin= vin,vmax=vax,cmap=plt.cm.get_cmap('nipy_spectral'))
plt.colorbar()
plt.plot(wh_p[::-1][:500],theta_res - 1 -(fl[wh_wl][:500]-np.amin(fl[wh_wl]))*25000,'k',linewidth=2)
plt.yticks(np.linspace(0,84,8))
plt.xticks(stick)

val = np.array([3,4,5,6,7,8,9])
stick = np.array([])
for i in range(val.size) :
    wh, = np.where(wl[obs_range][::-1] >= val[i])
    stick = np.append(stick,np.array([wh[0]]))

stick = stick - 500

plt.figure(4,figsize=(14,8))
plt.imshow(np.transpose(alpha[0,::-1,:][500:,:]) - np.transpose(alpha[1,::-1,:][500:,:]), aspect ='auto',vmin= vin,vmax=vax,cmap=plt.cm.get_cmap('nipy_spectral'))
plt.colorbar()
plt.plot(wh_p[::-1][498:]-502,95-(fl[wh_wl][498:]-np.amin(fl[wh_wl]))*25000,'k',linewidth=2)
plt.yticks(np.linspace(0,84,8))
plt.xticks(stick)

'''

val = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9,3,4,5,6,7,8,9])
stick = np.array([])
for i in range(val.size) :
    wh, = np.where(wl[obs_range][::-1] >= val[i])
    stick = np.append(stick,np.array([wh[0]]))

plt.figure(3,figsize=(14,8))
plt.imshow(np.transpose(alpha[0,::-1,:]) - np.transpose(alpha[1,::-1,:]), origin='lower',aspect ='auto',vmin= vin,vmax=vax,cmap=plt.cm.get_cmap('nipy_spectral'))
plt.plot(wh_p[::-1],(fl[wh_wl]-np.amin(fl[wh_wl]))*fact,'k',linewidth=2)
plt.colorbar()
plt.yticks(np.linspace(0,theta_res,9))
plt.xticks(stick)
plt.ylim(0,theta_res-1)
plt.xlim(0,s[1])


plt.show()