import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
from pyfunction import ratio
sys.path.append('/Users/caldas/Pytmosph3R/ComnPyRouts/')
from pytransfert import atmospectre
from pyconstant import *

#file = np.array(['tot','nude','cont','sca','0.50_0.50_ZnS','0.50_0.50_KCl'])
#file = np.linspace(-0.06,0.06,21)
file = np.array(['6.3','6.3','1.1','1.2','1.3','1.4','1.5','6.3','6.3','6.3'])
col = np.array(['k','k','r','b','g','brown','c','y','purple','orange'])
lab = np.array(['Absorption totale','Absorption sans nuages','H2O','CH4','NH3','CO','CO2','Continuum','Rayleigh','Diffusion de Mie'])
lw = np.array([3,3,2,2,2,2,2,2,2,2])
ls = np.array(['-','--','-','-','-','-','-','-','-','-','-'])
path = '/Users/caldas/Pytmosph3R'
source = 'Source'
name_source = 'bin10'
Rp = 0.246384689105*R_J
Rs = 0.206470165349*R_S

x = np.array([-1.392,-1.557,-8.129,-3.488,-8.557,-1.948])
x2 = np.array([-1.125,-1.425,-2.245,-2.383,-2.102,-1.949])

n_species = np.array(['H2','He','H2O','CH4','N2','NH3','CO','CO2'])

M_species, M, x_ratio_species = ratio(n_species,10**x,IsoComp=True)
print M
M_species, M, x_ratio_species = ratio(n_species,10**x2,IsoComp=True)
print M

weight = np.zeros(file.size-1)

#Rs = 0.114*R_S
#Rp = 1.12101710877*R_T

print (Rp**2)/(Rs**2)

plt.grid(True)
plt.semilogx()
axes = plt.gca()
for i_f in range(file.size) :

    i_f = file.size-1-i_f

    '''
    if file[i_f] == 'tot' :
        I1 = np.load('%s/I/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_52125_contscatot_IR.npy'%(path))
        I2 = np.load('%s/I/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_52125_0.50_0.50_ZnS_IR.npy'%(path))
        I3 = np.load('%s/I/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_52125_0.50_0.50_KCl_IR.npy'%(path))
        I = I1*I2*I3
    else :
        I = np.load('%s/I/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_52125_%s_IR.npy'%(path,file[i_f]))
    '''
    #I = np.load('/Users/caldas/Desktop/Test_rot/I_dis_6.3_Trappist_3_64x48_0_3000_200_1627_%s_IR.npy'%(file[i_f]))
    if i_f != 0 and i_f != 1 and i_f != 7 and i_f != 8 and i_f != 9 :
        I = np.load('/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_nude_IR.npy'%(file[i_f]))
    else :
        if i_f == 0 :
            #I = np.load(
            #'/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy' % (
            #file[i_f]))
            I = np.load(
            '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_0.50_0.50_tot_IR.npy' % (
            file[i_f]))
        if i_f == 1 :
            #I = np.load(
            #    '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_nude_IR.npy' % (
            #        file[i_f]))
            I = np.load(
            '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy' % (
            file[i_f]))
        if i_f == 7 :
            I = np.load(
                '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_cont_IR.npy' % (
                    file[i_f]))
        if i_f == 8 :
            I = np.load(
                '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_sca_IR.npy' % (
                    file[i_f]))
        if i_f == 9 :
            I = np.load(
                '/Users/caldas/Desktop/Heterogeneities/GJ1214b/I_GJ1214b/I_dis_%s_GJ1214b_3_64x48_0_3000_100_69615_0.50_0.50_cloud_IR.npy' % (
                    file[i_f]))


    bande_sample = np.load('%s/%s/bande_sample_%s.npy'%(path,source,name_source))
    wh_b, = np.where((10000/bande_sample >= 0.6)*(10000/bande_sample <= 20.))

    R_eff_bar,R_eff,ratio_bar,ratR_bar,bande_bar,flux_bar,flux = atmospectre(I,bande_sample,Rs,Rp,69615.,0,\
                                                                                    False,False,True)

    plt.plot(1/bande_sample*10000,flux,linestyle='%s'%(ls[i_f]), color = col[i_f], linewidth=lw[i_f], label=lab[i_f])

    if i_f == 0 :
        flux_ref = flux
    if i_f == 3 :
        flux_ref_min = flux

print 1/weight*100/(np.nansum(1/weight))


plt.axis([0.6,20.,np.amin(flux_ref_min),np.amax(flux_ref)+0.0002])
plt.xticks(np.array([0.7, 0.8, 0.9, 1., 2., 3., 3., 4., 5., 6., 7.,8.,9.,10.,20.]))
plt.ylabel('Relative flux')
plt.xlabel('Wavelength (micron)')
plt.legend(loc=4)
plt.show()