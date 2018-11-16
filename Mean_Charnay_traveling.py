import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
from pyconstant import *
from pyfunction import *
from matplotlib import colors as mcolors

path = '/Users/caldas/Desktop/Heterogeneities/'

col = np.array(['k','g','r','b','m'])
#Rp = 0.246384689105
Rp = 1.7e+7
print 1700000./R_J
print Rp/R_J
h = 6961500.
r_step = 12000000./100.
theta_step = np.pi/48.

n_cross = np.array([2,3,5,6,7])
n_data = np.array([6,7,9,10,11])
n_species = np.array(['H2','He','H2O','CH4','N2','NH3','CO','CO2'])
n_spe_array = np.array([7,2,3,5,6,7])
i_spe_array = np.array([4,0,1,2,3,4])
n_spe_n_array = np.array(['T', 'H2O', 'CH4', 'NH3', 'CO', 'CO2'])

bande_sample = np.load('/Users/caldas/Pytmosph3R/Source/bande_sample_bin10.npy')

bande_sample = 10**4/bande_sample

obs_range, = np.where((bande_sample > 1)*(bande_sample < 10.))
print bande_sample[obs_range]
crossection = np.load('/Users/caldas/Pytmosph3R/Source/crossection_bin10.npy')
print np.shape(crossection)
T_sample = np.load('/Users/caldas/Pytmosph3R/Source/T_sample_bin10.npy')
P_sample = np.load('/Users/caldas/Pytmosph3R/Source/P_sample_bin10.npy')
I = np.load('I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy')

stud = 'Para_GJ1214b/Parameters/GJ1214b_data_convert_6961x64x48_lo0.00_0.5mu_0.5mu'

order_grid = np.load('Para_GJ1214b/Stitch/order_grid_96_64x48x6961_69615_0.00_0.00.npy')
dx_grid = np.load('Para_GJ1214b/Stitch/dx_grid_opt_96_64x48x6961_69615_0.00_0.00.npy')
s = np.shape(dx_grid)

print s
dat = np.load('%s.npy' % (stud))
print np.shape(dat)

dx_op = np.zeros((4,96,2 * s[2]))
T_op = np.zeros((4,96,s[2]))
P_op = np.zeros((4,96,s[2]))
x_op = np.zeros((4,n_cross.size,96,s[2]))

z = np.array([42,43,44,45,46,47])

I = np.load('I_Charnay.npy')

plt.figure(1)
plt.grid(True)
axes = plt.gca()
i_thet = 24

sh = np.shape(I)

points = np.zeros(sh[1])
points_real = np.zeros(sh[1])
points_real_int = np.zeros(sh[1])

for i_b in range(sh[0]) :

    #wh, = np.where((I[i_b,:,0] < 0.99) * (I[i_b,:,0] > 0.01))
    for i_thet in range(96):
        wh, = np.where(np.abs(I[i_b, :-1, i_thet] - I[i_b, 1:, i_thet]) > 0.005)

        #points[wh] += 1/96.
        points[wh] += (I[i_b, wh + 1, i_thet] + I[i_b, wh - 1, i_thet]) / (2 * 96.)

        wh_b, = np.where(i_b == obs_range)

        if wh_b.size != 0:
            #points_real[wh] += 1/96.
            points_real[wh] += (I[i_b,wh+1,i_thet] + I[i_b,wh-1,i_thet])/(2 * 96.)

for i_r in range(sh[1]) :
    points_real_int[i_r] = np.nansum(points_real[:i_r])

np.save('Distrib_Limb.npy',points_real/obs_range.size)

plt.figure(1)
plt.plot(points/sh[0],'^',markersize=8,color=col[0])
plt.plot(points_real/obs_range.size,'o',markersize=8,color=col[0])
plt.show()


for i_z in range(4) :

    T_mean = (np.mean(dat[1,0,z[i_z],24,16]) + np.mean(dat[1,0,z[i_z],24,48]))/2.
    print T_mean

for i_z in range(4) :

    wh = np.array([z[i_z]])

    for i_thet in range(96) :

        range_op, = np.where(dx_grid[wh[0], i_thet, :] > 0)

        data = np.zeros(np.shape(dat))

        for i_lat in range(49):
            for i_long in range(65):
                if i_thet <= 48:
                    lat = i_thet
                    long = 16
                else:
                    lat = 96 - i_thet
                    long = 48
                data[:, :, :, i_lat, i_long] = dat[:, :, :, lat, long]

        for i_l in range(range_op.size):
            if i_l != 0:
                dx_op[i_z,i_thet, 2 * i_l] = np.nansum(dx_grid[wh[0], 24, 0:i_l])
            dx_op[i_z,i_thet, 2 * i_l + 1] = np.nansum(dx_grid[wh[0], 24, 0:i_l + 1])
            T_op[i_z,i_thet, i_l] = data[1, 0, order_grid[0, wh[0], 24, i_l], order_grid[1, wh[0], 24, i_l], order_grid[
                2, wh[0], 24, i_l]]
            P_op[i_z,i_thet, i_l] = data[0, 0, order_grid[0, wh[0], 24, i_l], order_grid[1, wh[0], 24, i_l], order_grid[
                2, wh[0], 24, i_l]]
            x_op[i_z,:,i_thet, i_l] = data[
                n_data, 0, order_grid[0, wh[0], 24, i_l], order_grid[1, wh[0], 24, i_l], order_grid[
                    2, wh[0], 24, i_l]]

        dx_op[i_z,i_thet] -= np.nansum(dx_grid[wh[0], i_thet, 0:range_op.size]) / 2.

del data, order_grid

X_MEAN = np.zeros((4,n_spe_array.size,96,obs_range.size))
TOT = np.zeros((4,n_spe_array.size,96,obs_range.size))

for i_z in range(4) :

    wh = np.array([z[i_z]])

    tau_op = np.zeros((96, s[2], obs_range.size))
    tau_op_mol = np.zeros((n_cross.size, 96, s[2], obs_range.size))

    for i_thet in range(96):

        range_op, = np.where(dx_grid[wh[0], i_thet, :] > 0)
        range_max = range_op[range_op.size - 1]+1

        res, c_grid, i_grid = interp2olation(np.log10(P_op[i_z,i_thet,range_op]) - 2, T_op[i_z,i_thet,range_op], P_sample, T_sample,
                                                 crossection[0,:, :, 0])

        c1, c2, c3, c4 = c_grid[:,0],c_grid[:,1],c_grid[:,2],c_grid[:,3]
        i1, i2, i3, i4 = i_grid[:,0],i_grid[:,1],i_grid[:,2],i_grid[:,3]

        del res, c_grid, i_grid

        for i_s in range(n_cross.size):
            cross = crossection[i_s, :, :, :]
            m = 0

            for i_b in obs_range :
                res = c3 * (c1 * cross[i2, i4, i_b] + c2 * cross[i1, i4, i_b]) + \
                      c2 * (c1 * cross[i2, i3, i_b] + c2 * cross[i1, i3, i_b])
                tau_op_mol[i_s, i_thet, :range_max, m] = res * P_op[i_z,i_thet,range_op] / (R_gp * T_op[i_z,i_thet,range_op]) * N_A * x_op[i_z,i_s,i_thet,range_op] * dx_grid[wh[0],i_thet,range_op]
                tau_op[i_thet, :range_max, m] += tau_op_mol[i_s, i_thet, :range_max, m]
                m += 1


    for i_t in range(n_spe_n_array.size):
        n_spe_n = n_spe_n_array[i_t]
        i_spe = i_spe_array[i_t]
        n_spe = n_spe_array[i_t]

        x_mean = np.zeros((96,obs_range.size))
        tot = np.zeros((96,obs_range.size))
        x_MEAN = 0

        for i_thet in range(96) :

            range_op, = np.where(dx_grid[wh[0], i_thet, :] > 0)
            range_max = range_op[range_op.size - 1] + 1
            m = 0

            for i_b in obs_range :
                tau_op_tot = np.nansum(tau_op[i_thet,:range_max,m])
                tau_op_mol_tot = np.nansum(tau_op_mol[i_spe,i_thet,:range_max,m])
                tau_pond = np.nansum(tau_op[i_thet,:range_max,m]*tau_op_mol[i_spe,i_thet,:range_max,m])
                if n_spe_n != 'T' :
                    x_mean[i_thet,m] = np.nansum(x_op[i_z,i_spe, i_thet,range_op]*tau_op[i_thet,:range_max,m]*tau_op_mol[i_spe,i_thet,:range_max,m])/(tau_pond)
                    if tau_op_mol_tot == 0. or tau_op_tot == 0 :
                        x_mean[i_thet,m] = 0.
                        tot[i_thet, m] = 0.
                    else :
                        tot[i_thet,m] = tau_op_mol_tot/tau_op_tot
                else :
                    x_mean[i_thet, m] = np.nansum(T_op[i_z,i_thet,range_op] * tau_op[i_thet, :range_max, m] / tau_op_tot)
                    tot[i_thet,m] = 1

                m += 1

        tot_tot = np.nansum(tot)

        x_MEAN = np.nansum(x_mean * tot) / tot_tot

        if n_spe_n != 'T' :
            print n_spe_n_array[i_t], np.log10(x_MEAN)
        else :
            print n_spe_n_array[i_t], x_MEAN

        X_MEAN[i_z,i_t] = x_mean
        TOT[i_z,i_t] = tot