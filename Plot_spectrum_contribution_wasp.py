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
file = np.array(['6.3','6.3','6.3','6.3'])
col = np.array(['k','k','r','b'])
lab = np.array(['Absorption totale','Absorption moleculaire','Continuum','Rayleigh'])
lw = np.array([3,3,2,2])
ls = np.array(['-','--','-','-'])
stus = np.array(['contscatot','nude','cont','sca'])
path = '/Users/caldas/Pytmosph3R'
source = 'Source'
name_source = 'bin10'
Rs = 1.458 * R_S
Rp = 1.865*R_J

weight = np.zeros(file.size-1)

#Rs = 0.114*R_S
#Rp = 1.12101710877*R_T

print (Rp**2)/(Rs**2)

plt.grid(True)
plt.semilogx()
axes = plt.gca()
for i_f in range(file.size) :

    if lab[i_f] != 'Absorption totale' :
        I = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518213_%s_IR.npy'%(stus[i_f]))

    else :
        I = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518213_la0.000_lo0.000_%s_IR.npy' % (stus[i_f]))


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