import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
from pyconstant import *
from pyfunction import *
from matplotlib import colors as mcolors

Already_saved = False
Already_tau = False

div = 100
#Rs = 0.206470165349*R_S
Rs = 1.458 * R_S
#Rp = 1.7e+7
Rp = 1.865*R_J

exo = 'WASP121'
name = 'H2O_CO2_var_2'

#n_cross = np.array([2,3,5,6,7])
#n_data = np.array([6,7,9,10,11])
#n_species = np.array(['H2','He','H2O','CH4','N2','NH3','CO','CO2'])
#n_species_active = np.array(['H2O','CH4','NH3','CO','CO2'])

n_cross =  np.array([0,3])
n_data = np.array([4,5])
n_species = np.array(['H2','He','H2O','CO'])
n_species_active = np.array(['H2O','CH4','NH3','CO','CO2'])

bande_sample = np.load('/Users/caldas/Pytmosph3R/Source/bande_sample_bin10.npy')
wl = 10**4/bande_sample

obs_range, = np.where((wl > 1.)*(wl < 10.))
crossection = np.load('/Users/caldas/Pytmosph3R/Source/crossection_bin10.npy')
crossection = crossection[:,:,:,obs_range]*1.e-4
class star:
    def __init__(self):
        #self.radius = 0.206470165349*R_S
        #self.temperature = 3000.
        #self.distance = 42.4*9.461e+15
        self.radius = 1.458 * R_S
        self.temperature = 6460.
        self.distance = 85.*9.461e+15
detection = JWST()
resolution = ''
int_lambda = np.zeros((2,bande_sample.size))
for i_bande in range(bande_sample.size) :
    if i_bande == 0 :
        int_lambda[0,i_bande] = bande_sample[0]
        int_lambda[1,i_bande] = (bande_sample[i_bande+1]+bande_sample[i_bande])/2.
    elif i_bande == bande_sample.size - 1 :
        int_lambda[0,i_bande] = (bande_sample[i_bande-1]+bande_sample[i_bande])/2.
        int_lambda[1,i_bande] = bande_sample[bande_sample.size-1]
    else :
        int_lambda[0,i_bande] = (bande_sample[i_bande-1]+bande_sample[i_bande])/2.
        int_lambda[1,i_bande] = (bande_sample[i_bande+1]+bande_sample[i_bande])/2.
int_lambda = np.sort(10000./int_lambda[::-1])
noise = stellar_noise(star(),detection,int_lambda,resolution,Total=False)
noise = noise[::-1]
bande_sample = bande_sample[obs_range]
noise = noise[obs_range]

print noise

sh_cr = np.shape(crossection)
print sh_cr

T_sample = np.load('/Users/caldas/Pytmosph3R/Source/T_sample_bin10.npy')
P_sample = np.load('/Users/caldas/Pytmosph3R/Source/P_sample_bin10.npy')
#I = np.load('I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy')
#I = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_516316_nude_IR.npy')
I = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518212_nude_IR.npy')

I = I[obs_range,:,:]
#r_step = 69615.
#r_step = 516316.
r_step = 518213.
sh = np.shape(I)
print sh

#stud = 'Para_GJ1214b/Parameters/GJ1214b_data_convert_6961x64x48_lo0.00_0.5mu_0.5mu'
#dx_grid = np.load('Para_GJ1214b/Stitch/dx_grid_opt_96_64x48x6961_69615_0.00_0.00.npy')
#stud = 'Para_WASP/Parameters/WASP121b_40_data_convert_51631x128x64_lo0.00_0.5mu'
#dx_grid = np.load('Para_WASP/Stitch/dx_grid_opt_128_128x64x51631_516316_0.00_0.00.npy')
stud = 'Para_WASP/Parameters/WASP121b_40_data_convert_518212x128x64_lo0.00_0.5mu'
dx_grid = np.load('Para_WASP/Stitch/dx_grid_opt_128_128x64x51821_518212_0.00_0.00.npy')
dx_grid = dx_grid[:sh[1]]
s = np.shape(dx_grid)
print s

indice = np.zeros((s[1],sh[0]),dtype=np.int)
indice_board = np.zeros((2,s[1],sh[0]),dtype=np.int)
val_board = np.zeros((2,s[1],sh[0]))

for i_b in range(sh[0]) :

    for i_thet in range(s[1]) :

        wh_I, = np.where(I[i_b, :, i_thet] >= np.exp(-1))
        if wh_I.size != 0 :
            indice[i_thet,i_b] = wh_I[0] - 1
            if I[i_b,wh_I[0]-1,i_thet] < 0.1 :
                indice[i_thet, i_b] = wh_I[0]
        else :
            print 'There is a problem with bande ', wl[i_b]
        wh_95, = np.where(I[i_b, :, i_thet] >= 0.95)
        wh_05, = np.where(I[i_b, :, i_thet] >= 0.05)
        if wh_95.size != 0 :
            indice_board[1,i_thet,i_b] = wh_95[0]
            val_board[1,i_thet,i_b] = I[i_b,wh_95[0],i_thet]
        if wh_05.size != 0 :
            indice_board[0,i_thet,i_b] = wh_05[0]
            val_board[0,i_thet,i_b] = I[i_b,wh_05[0],i_thet]

ind_min = np.int(np.amin(indice))
ind_max = np.int(np.amax(indice))+1
print ind_min, ind_max

range_m = np.zeros((ind_max-ind_min,s[1]),dtype=np.int)

for i_thet in range(s[1]):
    already = np.array([])
    for i_b in range(sh[0]):

        i_z = indice[i_thet, i_b]

        wh_iz, = np.where(already == i_z)

        if wh_iz.size == 0:
            already = np.append(already, np.array([i_z]))
            rangee, = np.where(dx_grid[i_z, i_thet, :] > 0)
            range_m[i_z-ind_min, i_thet] = rangee.size

np.save('Files_penetration_duo/%s/indice_%s.npy'%(exo,name),indice)
np.save('Files_penetration_duo/%s/indice_board_%s.npy'%(exo,name),indice_board)
np.save('Files_penetration_duo/%s/val_board_%s.npy'%(exo,name),val_board)

if Already_saved == False :

    data = np.load('%s.npy' % (stud))

    #order_grid = np.load('Para_GJ1214b/Stitch/order_grid_96_64x48x6961_69615_0.00_0.00.npy')
    #order_grid = np.load('Para_WASP/Stitch/order_grid_128_128x64x51631_516316_0.00_0.00.npy')
    order_grid = np.load('Para_WASP/Stitch/order_grid_128_128x64x51821_518212_0.00_0.00.npy')

    order_grid = order_grid[:sh[1]]

    T_op = np.zeros((sh[0], s[1], s[2]))
    P_op = np.zeros((sh[0], s[1], s[2]))
    x_op = np.zeros((sh[0], n_species_active.size, s[1], s[2]))

    bar = ProgressBar(s[1], 'Progression/Parametrisation : ')

    for i_thet in range(s[1]) :

        for i_b in range(sh[0]) :

            i_zt = indice[i_thet, i_b]
            i_z = i_zt - ind_min
            range_max = range_m[i_z,i_thet]

            T_op[i_b,i_thet, :range_max] = data[1, 0, order_grid[0, i_zt, i_thet, :range_max], order_grid[1, i_zt, i_thet, :range_max], order_grid[
                    2, i_zt, i_thet, :range_max]]
            P_op[i_b,i_thet, :range_max] = data[0, 0, order_grid[0, i_zt, i_thet, :range_max], order_grid[1, i_zt, i_thet, :range_max], order_grid[
                    2, i_zt, i_thet, :range_max]]
            for i_s in range(n_data.size) :
                x_op[i_b,i_s,i_thet, :range_max] = data[n_data[i_s], 0, order_grid[0, i_zt, i_thet, :range_max], order_grid[1, i_zt, i_thet, :range_max], order_grid[
                    2, i_zt, i_thet, :range_max]]
        bar.animate(i_thet+1)

    np.save('Files_penetration_duo/%s/T_op_%s.npy' % (exo, name), T_op)
    np.save('Files_penetration_duo/%s/P_op_%s.npy' % (exo, name), P_op)
    np.save('Files_penetration_duo/%s/x_op_%s.npy' % (exo, name), x_op)

    del data, order_grid

else:

    T_op = np.load('Files_penetration_duo/%s/T_op_%s.npy' % (exo, name))
    P_op = np.load('Files_penetration_duo/%s/P_op_%s.npy' % (exo, name))
    x_op = np.load('Files_penetration_duo/%s/x_op_%s.npy' % (exo, name))

print 'Coefficient calculation is started'

if Already_tau == False :
    c31 = np.zeros((ind_max-ind_min,s[1],s[2]))
    c32 = np.zeros((ind_max-ind_min,s[1],s[2]))
    c41 = np.zeros((ind_max-ind_min,s[1],s[2]))
    c42 = np.zeros((ind_max-ind_min,s[1],s[2]))
    i1 = np.zeros((ind_max-ind_min,s[1],s[2]),dtype=np.int)
    i2 = np.zeros((ind_max-ind_min,s[1],s[2]),dtype=np.int)
    i3 = np.zeros((ind_max-ind_min,s[1],s[2]),dtype=np.int)
    i4 = np.zeros((ind_max-ind_min,s[1],s[2]),dtype=np.int)
    fac = np.zeros((s[1],sh[0],s[2]))

    bar = ProgressBar(s[1], 'Progression/Coefficients calculations : ')

    for i_thet in range(s[1]) :
        already = np.array([])
        for i_b in range(sh[0]) :

            i_zt = indice[i_thet, i_b]
            i_z = i_zt - ind_min

            range_max = range_m[i_z, i_thet]

            wh_z, = np.where(already == i_zt)

            if wh_z.size == 0 :

                res, c_grid, i_grid = interp2olation(np.log10(P_op[i_b, i_thet, :range_max]) - 2, T_op[i_b, i_thet, :range_max], P_sample, T_sample, crossection[0,:, :, 0])

                c31[i_z, i_thet, :range_max] = c_grid[:, 0]*c_grid[:, 2]
                c32[i_z, i_thet, :range_max] = c_grid[:, 1]*c_grid[:, 2]
                c41[i_z, i_thet, :range_max] = c_grid[:, 0]*c_grid[:, 3]
                c42[i_z, i_thet, :range_max] = c_grid[:, 1]*c_grid[:, 3]
                i1[i_z, i_thet, :range_max] = i_grid[:, 0]
                i2[i_z, i_thet, :range_max] = i_grid[:, 1]
                i3[i_z, i_thet, :range_max] = i_grid[:, 2]
                i4[i_z, i_thet, :range_max] = i_grid[:, 3]

                already = np.append(already, np.array([i_zt]))

            fac[i_thet,i_b,:range_max] = P_op[i_b, i_thet, :range_max] / (R_gp * T_op[i_b, i_thet, :range_max]) * N_A * dx_grid[i_zt,i_thet,:range_max]

        bar.animate(i_thet+1)

    tau_op_mol_save = np.zeros((n_species_active.size, s[1], s[2], sh[0]))

    bar = ProgressBar(obs_range.size, 'Progression/Optical depth calculations : ')

    for i_b in range(sh[0]):

        for i_thet in range(s[1]):

            i_zt = indice[i_thet, i_b]
            i_z = i_zt - ind_min

            range_max = range_m[i_z,i_thet]

            ia = i1[i_z,i_thet,:range_max]
            ib = i2[i_z,i_thet,:range_max]
            ic = i3[i_z,i_thet,:range_max]
            id = i4[i_z,i_thet,:range_max]

            cross31 = crossection[:, ib, id, i_b]
            cross32 = crossection[:, ia, id, i_b]
            cross41 = crossection[:, ib, ic, i_b]
            cross42 = crossection[:, ia, ic, i_b]

            fac_i = fac[i_thet,i_b,:range_max]

            m = 0
            for i_s in n_cross :

                res = c31[i_z,i_thet,:range_max] * cross31[i_s] \
                    + c32[i_z,i_thet,:range_max] * cross32[i_s] \
                    + c41[i_z,i_thet,:range_max] * cross41[i_s] \
                    + c42[i_z,i_thet,:range_max] * cross42[i_s]
                tau_op_mol_save[m, i_thet, :range_max, i_b] = res * x_op[i_b, m, i_thet, :range_max] * fac_i
                m += 1
        bar.animate(i_b + 1)

    np.save('Files_penetration_duo/%s/tau_op_mol_save_%s.npy' % (exo, name), tau_op_mol_save)

else:

    tau_op_mol_save = np.load('Files_penetration_duo/%s/tau_op_mol_save_%s.npy' % (exo, name))

alpha = np.zeros((2, sh[0], s[1]))

print 'Tau calculation is started'

correction_b = np.array([])
correction_t = np.array([])
correction = np.array([])

bar = ProgressBar(obs_range.size, 'Progression/Penetration angle calculations : ')

for i_b in range(sh[0]) :

    for i_thet in range(s[1]):

        i_zt = indice[i_thet,i_b]
        i_z = i_zt - ind_min
        range_max = range_m[i_z, i_thet]

        D = np.nansum(dx_grid[i_zt, i_thet, :]) / 2.
        dx_cal = dx_grid[i_zt,i_thet,:range_max]

        delta_R = Rs ** 2 * noise[i_b] / (2 * (Rp + (i_zt + 0.5) * r_step))
        a = (val_board[1, i_thet, i_b] - val_board[0, i_thet, i_b]) / (
                    (indice_board[1, i_thet, i_b] - indice_board[0, i_thet, i_b]) * r_step)
        delta_I = a * delta_R

        tau_op_tot = np.nansum(tau_op_mol_save[:, i_thet, :range_max, i_b])

        if np.exp(-tau_op_tot) < 0.05 :
            print 'prob for ', i_thet, i_b
            print np.exp(-tau_op_tot)
            correction_b = np.append(correction_b,np.array([i_b]))
            correction_t = np.append(correction_b,np.array([i_b]))
            correction = np.append(correction, np.array(['d']))
        if np.exp(-tau_op_tot) > 0.95 :
            print 'prob for ', i_thet, i_b
            print np.exp(-tau_op_tot)
            correction_b = np.append(correction_b,np.array([i_b]))
            correction_t = np.append(correction_b,np.array([i_b]))
            correction = np.append(correction, np.array(['u']))

        t = 1
        while np.exp(-np.nansum(tau_op_mol_save[:,i_thet, t:range_max,i_b])) < np.exp(-tau_op_tot) + delta_I / 2.:
            t += 1

        u = 1
        while np.exp(-np.nansum(tau_op_mol_save[:,i_thet, :range_max-u+1,i_b])) < np.exp(-tau_op_tot) + delta_I / 2.:
            u += 1

        tau_star = np.nansum(tau_op_mol_save[:,i_thet,t-1:range_max,i_b])
        tau_obs = np.nansum(tau_op_mol_save[:,i_thet, :range_max-u+2,i_b])

        tau_op_mol_star = np.nansum(tau_op_mol_save[:,i_thet, t-1, i_b]) / np.float(div)
        tau_op_mol_obs = np.nansum(tau_op_mol_save[:, i_thet, range_max-u+1, i_b]) / np.float(div)

        m = 1
        while np.exp(-np.nansum(tau_star - m * tau_op_mol_star)) < np.exp(-tau_op_tot) + delta_I/2. :
            m += 1

        n = 1
        while np.exp(-np.nansum(tau_obs - n * tau_op_mol_obs)) < np.exp(-tau_op_tot) + delta_I/2. :
            n += 1

        dx_moins = D - np.nansum(dx_cal[:t-1]) - dx_cal[t-1] / np.float(div) * (m-1)
        dx_plus = D - np.nansum(dx_cal[range_max-u+2:]) - dx_cal[range_max-u+1] / np.float(div) * (n-1)

        alpha[0, i_b, i_thet] = np.arctan2(dx_moins, Rp + (i_zt + 0.5) * r_step)*360/(2*np.pi)
        alpha[1, i_b, i_thet] = np.arctan2(dx_plus, Rp + (i_zt + 0.5) * r_step)*360/(2*np.pi)

    bar.animate(i_b + 1)

np.save('Files_penetration_duo/%s/alpha_%s.npy'%(exo,name),alpha)
np.save('Files_penetration_duo/%s/correction_%s.npy'%(exo,name),correction)
np.save('Files_penetration_duo/%s/correction_t_%s.npy'%(exo,name),correction_t)
np.save('Files_penetration_duo/%s/correction_b_%s.npy'%(exo,name),correction_b)

plt.imshow(alpha[0,:,:], aspect ='auto')
plt.colorbar()
plt.show()

plt.imshow(alpha[0,:,:], aspect ='auto')
plt.colorbar()
plt.show()