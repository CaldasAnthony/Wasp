import numpy as np
from pyfunction import *
from pyconstant import *
from pydataread import *
from pytransfert import *


#Rp = 15.*R_T
#Rs = 1.155*R_S
#r_step = 69000.
#r_step = 108000.
r_step = 69615.
duo_To = '500_1000'
duo_T = '1000_1800'
#Rp = 1.7e+7
#Rs = 0.206470165349*R_S
#Ts = 3000.
#d_al = 42.*9.461e+15
Rp = 1.865*R_J
Rs = 1.458 * R_S
Ts = 6460.
d_al = 85.*9.461e+15

w_min = 1.
w_max = 10.

path = '/Users/caldas/Pytmosph3R/'
name_source = 'Source'
source = 'bin10'
resolution = ''
detection = JWST()

class star:
    def __init__(self):
        self.radius = Rs
        self.temperature = Ts
        self.distance = d_al

bande_sample = np.load("%s%s/bande_sample_%s.npy" % (path, name_source, source))
bande_sample = np.delete(bande_sample, [0])
#wl = 1 / bande_sample * 10000
#wh, = np.where((wl >= w_min)*(wl <=w_max))
#bande_sample = bande_sample[wh]
int_lambda = np.zeros((2, bande_sample.size))
bande_sample = np.sort(bande_sample)


#I = np.load('I_Charnay.npy')
#save_ad = 'GJ1214b_Limb'

#Ia = np.load('I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.npy')
Ia = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518212_nude_IR.npy')
#Ia = np.load('I_WASP/I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518213_la0.000_lo0.000_contscatot_IR.npy')
save_ad = 'WASP121_2'

#I = Ia[wh,:,:]
I = Ia
#save_ad = 'GJ1214b_Est_no'

#I = np.zeros((3000,100,80))
#I[:,:,0:64] = Ia[:,:,0:64]
#I[:,:,64:] = Ia[:,:,80:]
#save_ad = 'GJ1214b_Ouest'

#plt.imshow(Ia[wh[0],:,:],aspect='auto',origin = 'lower')
#plt.colorbar()
#plt.show()

#I = Ia[:,:,64:80]
#save_ad = 'GJ1214b_Est'

#plt.imshow(I[wh[0],:,:],aspect='auto',origin = 'lower')
#plt.colorbar()
#plt.show()

if resolution == '':
    int_lambda = np.zeros((2, bande_sample.size))
    for i_bande in range(bande_sample.size):
        if i_bande == 0:
            int_lambda[0, i_bande] = bande_sample[0]
            int_lambda[1, i_bande] = (bande_sample[i_bande + 1] + bande_sample[i_bande]) / 2.
        elif i_bande == bande_sample.size - 1:
            int_lambda[0, i_bande] = (bande_sample[i_bande - 1] + bande_sample[i_bande]) / 2.
            int_lambda[1, i_bande] = bande_sample[bande_sample.size - 1]
        else:
            int_lambda[0, i_bande] = (bande_sample[i_bande - 1] + bande_sample[i_bande]) / 2.
            int_lambda[1, i_bande] = (bande_sample[i_bande + 1] + bande_sample[i_bande]) / 2.
    int_lambda = np.sort(10000. / int_lambda[::-1])
else:
    int_lambda = np.sort(10000. / bande_sample[::-1])

noise = stellar_noise(star(), detection, int_lambda, resolution)
noise = noise[::-1]

bande_sample = np.load("%s%s/bande_sample_%s.npy" % (path, name_source, source))
#wh, = np.where((wl >= w_min)*(wl <=w_max))
#bande_sample = bande_sample[wh]

flux_script(path,name_source,source,save_ad,I,noise,Rs,Rp,r_step,bande_sample,Kcorr=False,Middle=True,Noise=False)

rec = atmospectre(I,bande_sample,Rs,Rp,r_step,0,False,False,True)

rec[6]

plt.figure(1)
plt.semilogx()

data = open('GJ1214b_no.dat', 'r')
I = data.readlines()

fl = np.zeros(len(I))
wl = np.zeros(len(I))
for i in range(len(I)):
    fl[i] = line_search(I[i])[1]
    wl[i] = line_search(I[i])[0]

plt.plot(wl,fl)

data = open('I_GJ1214b/I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR.dat', 'r')
I = data.readlines()

fl = np.zeros(len(I))
wl = np.zeros(len(I))
for i in range(len(I)):
    fl[i] = line_search(I[i])[1]
    wl[i] = line_search(I[i])[0]

plt.plot(1/bande_sample*10**4,rec[6])
plt.plot(wl,fl)
plt.show()
