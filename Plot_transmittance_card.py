from pyfunction import *

#I_tot = np.load('/Users/caldas/Pytmosph3R/I/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_52125_contscatot_IR.npy')
I_tot = np.load('/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy')
#I_tot = np.load('/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_0.50_0.50_tot_IR.npy')
sh = np.shape(I_tot)

n_layers = 80

I = I_tot[:,:n_layers,:]

"""
for i in range(3000) :
    plt.plot(I[i,:,0])
plt.show()

print sh

I = np.zeros(sh)
for i_theta in range(sh[2]) :
    I[:,:,i_theta] = I_tot[:,:,(sh[2]/2+i_theta)%(sh[2])]
"""

r_step = 69615.
h = n_layers*r_step

theta_number = 96
name = 'GJ1214b'
Rp = 0.246384689105*R_J
Rp = 1.7e+7
Rs = 0.206470165349*R_S

param = 20000
factor = 10.

#wl = 0.95
wl = 2.70
bande_sample = np.load('/Users/caldas/Pytmosph3R/Source/bande_sample_bin10.npy')

I_png = atmosphere_plot(I,h,param,factor,r_step,theta_number,wl,bande_sample,name,Rp,Rs,0,0,'polar',False,False,False)

np.save('I_NH3_2.70.npy',I_png)