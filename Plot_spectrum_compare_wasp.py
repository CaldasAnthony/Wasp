import numpy as np
import matplotlib.pyplot as plt

import sys

sys.path.append('/Users/caldas/Pytmosph3R/ComnPyRouts/')
#from netCDF4 import Dataset
from pyfunction import *
from pyconstant import *
from pytransfert import generator_1D
from pytransfert import atmospectre

cont_tot = np.array(['H2-He_2011.cia','H2-He_2011.cia'])
cont_species = np.array(['H2','He'])
cont_associations = np.array(['h2h2','h2he'])
class continuum :
    def __init__(self) :
        self.number = cont_tot.size
        self.associations = cont_associations
        self.species = cont_species
        self.file_name = cont_tot

directory = '/Users/caldas/Pytmosph3R/'
path = '/Users/caldas/Heterogeneities/WASP121b/'
#name_data = 'WASP121b_40_data_convert_51821x128x64_lo0.00_0.5mu'
name_data = 'WASP121b_40_data_convert_51631x128x64_lo0.00_0.5mu'
#name_data = 'GJ1214b_data_convert_52126448'
#name_I = 'I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_518213_la0.000_lo0.000_contscatot_IR'
name_I = 'I_dis_6.3_WASP121b_40_3_128x64_0_3000_100_516316_la0.000_lo0.000_contscatot_IR'
#name_I = 'I_dis_6.3_GJ1214b_3_64x48_0_3000_100_69615_la0.000_lo0.000_contscatot_IR'
name_file = 'GJ'
domain = 'HR'

T_comp = np.load('%sSource/T_comp_WASP121b_40.npy'%(directory))
P_comp = np.load('%sSource/P_comp_WASP121b_40.npy'%(directory))
x_species = np.load('%sSource/x_species_comp_WASP121b_40.npy'%(directory))
n_species = np.array(['H2','He','H2O','CO'])
n_species_cross = np.array(['H2O','CO'])
n_species_active = np.array(['H2O','CO'])
x_ratio_species_active = np.array([0.01,0.01])
M_species, M, x_ratio_species = ratio(n_species,x_ratio_species_active,IsoComp=False)
ind_cross, ind_active = index_active (n_species,n_species_cross,n_species_active)

Charnay = True
D1_extrapolate = False
D1 = False
D1_mean = False
Noise = True
No_extra = False
Cl = False
E_O = False

m_species = np.array([])
c_species = np.array([])
dim_bande = 3000
Rp = 1.865*R_J
Mp = 1.184*M_J
g0 = G*Mp/(Rp**2)
print g0

Rs = 1.458 * R_S
Ts = 6450.
d_al = 85.*9.461e+15
phi_rot = 0.00

#h = 51821300.
#h0 = 51821300.
h = 51631600.
h0 = 51631600.
r_step = h/100.
n_layers = 100

cross = np.load('%sSource/crossection_bin10.npy'%(directory))
cross = cross[ind_cross]
k_cont = continuum()
Q_cloud = np.load('%sSource/Q_KCl_ZnS_GJ1214b.npy'%(directory))
bande_cloud = np.load('%sSource/bande_cloud_GJ1214b.npy'%(directory))
r_cloud = np.load('%sSource/radius_cloud_GJ1214b.npy'%(directory))
rho_p = np.array([])
r_eff = np.array([])
P_sample = np.load('%sSource/P_sample_bin10.npy'%(directory))
T_sample = np.load('%sSource/T_sample_bin10.npy'%(directory))
Q_sample = np.array([])
bande_sample = np.load('%sSource/bande_sample_bin10.npy'%(directory))


reso_theta = 128
reso_lat = 64
reso_long = 128

#data_convert = np.load('%sFiles/Parameters/%s.npy'%(directory,name_data))
data_convert = np.load('%sPara_WASP/Parameters/%s.npy'%(path,name_data))
sh_d = np.shape(data_convert)
print sh_d

#file = Dataset("/Users/caldas/Bureautique/Simu/GJ1214b/diagfi.nc")
#variables = file.variables
#sh = np.shape(variables['p'][:])
#pmin = np.amin(variables['p'][5,:,:,:])
#pmax = np.amax(variables['p'][5,:,:,:])

if D1 == True :
    diagfi = np.load('GJ1214b.npy')
    diagfi_surf = np.load('GJ1214b_surf.npy')
    sh = np.shape(diagfi[0])
    pmin = np.amin(diagfi[0,5])
    pmax = np.amax(diagfi[0,5])
    print pmin

    alt_array = np.zeros((reso_theta,sh[1]+1))

flux_extra_dis = np.zeros((reso_theta,dim_bande))
flux_extra_int = np.zeros((reso_theta,dim_bande))
flux_dis = np.zeros((reso_theta,dim_bande))
flux_int = np.zeros((reso_theta,dim_bande))
flux_mean = np.zeros(dim_bande)
i_max = 0

I_Charnay = np.zeros((dim_bande,n_layers,reso_theta))

ratio_HeH2 = 0.0000001

if E_O == True :
    ran = np.array([24,72])
else :
    ran = range(reso_theta)

for i_theta in ran :
    if i_theta > reso_lat :
        i_lat = reso_theta - i_theta
        i_long = reso_long/4
    else :
        i_lat = i_theta
        i_long = 3*reso_long/4

    if D1 == True :

        data = np.zeros((sh_d[0],sh[1]+1))
        '''
        data[0,1:] = variables['p'][5,:,i_lat,i_long]
        data[1,1:] = variables['temp'][5,:,i_lat,i_long]
        data[2,1:] = variables['gen_cond'][5,:,i_lat,i_long]
        data[3,1:] = variables['gen_cond2'][5,:,i_lat,i_long]
        data[0,0] = variables['ps'][5,i_lat,i_long]
        data[1,0] = variables['tsurf'][5,i_lat,i_long]
        data[2,0] = variables['gen_cond_surf'][5,i_lat,i_long]
        data[3,0] = variables['gen_cond2_surf'][5,i_lat,i_long]
        '''
        data[0, 1:] = diagfi[0, 5, :, i_lat, i_long]
        data[1, 1:] = diagfi[1, 5, :, i_lat, i_long]
        data[2, 1:] = diagfi[2, 5, :, i_lat, i_long]
        data[3, 1:] = diagfi[3, 5, :, i_lat, i_long]
        data[0, 0] = diagfi_surf[0, 5, 0, i_lat, i_long]
        data[1, 0] = diagfi_surf[1, 5, 0, i_lat, i_long]
        data[2, 0] = diagfi_surf[2, 5, 0, i_lat, i_long]
        data[3, 0] = diagfi_surf[3, 5, 0, i_lat, i_long]

        res, c_grid, i_grid = interp2olation_multi(data[0],data[1],P_comp,T_comp,x_species)
        data[6:sh_d[0]-1] = res[2:]
        for i in range(sh[1]+1) :
            data[4,i] = (1. - np.nansum(data[6:sh_d[0]-1,i]))/(ratio_HeH2 + 1.)
            data[5,i] = data[4,i]*ratio_HeH2

        data[sh_d[0]-1] = np.dot(M_species,data[4:sh_d[0]-1,:])

########################################################################################################################

    if D1 == True :

        alt_array[i_theta,0] = 0.
        for i_alt in range(1,sh[1]+1) :
            '''
            g = g0*1/(1+alt_array[i_theta,i_alt-1]/Rp)**2
            dz = R_gp*data[1,i_alt-1]/(data[sh_d[0]-1,i_alt-1]*g)*np.log(data[0,i_alt-1]/data[0,i_alt])
            '''
            a_z = -(1+alt_array[i_theta,i_alt-1]/Rp)*R_gp*(data[1,i_alt]-data[1,i_alt-1])/((data[sh_d[0]-1,i_alt]+data[sh_d[0]-1,i_alt-1])/2.*g0*\
                    np.log(data[1,i_alt]/data[1,i_alt-1]))*np.log(data[0,i_alt]/data[0,i_alt-1])
            dz = a_z*(1+alt_array[i_theta,i_alt-1]/Rp)/(1-a_z/Rp)
            alt_array[i_theta,i_alt] = alt_array[i_theta,i_alt-1]+dz

        h = alt_array[i_theta,sh[1]]
        print h

        I_dis, I_int = generator_1D(data,n_species,m_species,c_species,dim_bande,Rp,h,alt_array[i_theta],g0,cross,ind_active,k_cont,\
                     Q_cloud,P_sample,T_sample,Q_sample,bande_sample,bande_cloud,r_eff,r_cloud,rho_p,name_file,5,phi_rot,\
                     domain,ratio,directory,Tracer=False,Continuum=True,Scattering=True,Clouds=Cl,Kcorr=False,\
                     Optimal=False,Discret=True,Integral=True,Gravity=False,Molecular=False,Pressure=True,Script=False,Middle=True)

        I = np.zeros((dim_bande,sh[1]+1))
        I[:,0] = I_dis[:,0]
        I[:,sh[1]] = I_dis[:,sh[1]-2]
        I[:,1:sh[1]] = I_dis
        for i_b in range(dim_bande) :
            alpha = 2*integrate.trapz((Rp+alt_array[i_theta])*(1-I[i_b,:]),alt_array[i_theta])
            flux_dis[i_theta,i_b] = (Rp**2+alpha)/Rs**2
        I = np.zeros((dim_bande,sh[1]+1))
        I[:,0] = I_int[:,0]
        I[:,sh[1]] = I_int[:,sh[1]-2]
        I[:,1:sh[1]] = I_int
        for i_b in range(dim_bande) :
            alpha = 2*integrate.trapz((Rp+alt_array[i_theta])*(1-I[i_b,:]),alt_array[i_theta])
            flux_int[i_theta,i_b] = (Rp**2+alpha)/Rs**2

    #####

    pmin = np.amin(data_convert[:,0,:,:,:])
    da = data_convert[:,0,:,i_lat,i_long]
    wh, = np.where(da[0] >= pmin)
    da = da[:,wh]
    h = (wh.size-2)*r_step+r_step/2.
    siz = wh.size
    alt = np.zeros(wh.size)
    alt[1:wh.size] = np.arange(r_step/2.,h+r_step,r_step)
    alt[0] = 0
    if i_max < wh.size :
        i_max = wh.size
    print h

########################################################################################################################

    if Charnay == True :

        I_dis, I_int = generator_1D(da,n_species,m_species,c_species,dim_bande,Rp,h,r_step,g0,cross,ind_active,k_cont,\
                     Q_cloud,P_sample,T_sample,Q_sample,bande_sample,bande_cloud,r_eff,r_cloud,rho_p,name_file,5,phi_rot,\
                     domain,ratio,directory,Tracer=False,Continuum=True,Scattering=True,Clouds=Cl,Kcorr=False,\
                     Optimal=False,Discret=True,Integral=True,Gravity=False,Molecular=True,Pressure=False,Script=False,Middle=True)
        '''
        I_diss = np.zeros((np.shape(I_dis)[0],np.shape(I_dis)[1],1))
        I_diss[:,:,0] = I_dis
        R_eff_bar,R_eff,ratio_bar,ratR_bar,bande_bar,flux_bar,flux = atmospectre(I_diss,bande_sample,Rs,Rp,r_step,0,\
                                                                                    False,False,True)
        flux_extra_dis[i_theta,:] = flux
        I_intt = np.zeros((np.shape(I_int)[0],np.shape(I_int)[1],1))
        I_intt[:,:,0] = I_int
        R_eff_bar,R_eff,ratio_bar,ratR_bar,bande_bar,flux_bar,flux = atmospectre(I_intt,bande_sample,Rs,Rp,r_step,0,\
                                                                                    False,False,True)
        flux_extra_int[i_theta,:] = flux
        '''
        I_Charnay[:,:,i_theta] = I_dis

        I = np.zeros((dim_bande,wh.size))
        I[:,0] = I_dis[:,0]
        I[:,wh.size-1] = I_dis[:,wh.size-3]
        I[:,1:wh.size-1] = I_dis
        for i_b in range(dim_bande) :
            alpha = 2*integrate.trapz((Rp+alt)*(1-I[i_b,:]),alt)
            flux_extra_dis[i_theta,i_b] = (Rp**2+alpha)/Rs**2
        I = np.zeros((dim_bande,wh.size))
        I[:,0] = I_int[:,0]
        I[:,wh.size-1] = I_int[:,wh.size-3]
        I[:,1:wh.size-1] = I_int
        for i_b in range(dim_bande) :
            alpha = 2*integrate.trapz((Rp+alt)*(1-I[i_b,:]),alt)
            flux_extra_int[i_theta,i_b] = (Rp**2+alpha)/Rs**2

    print i_theta

    '''
    plt.semilogx()
    plt.grid(True)
    axes = plt.gca()
    if i_theta == 24 :
        plt.plot(1 / (100. * bande_sample) * 10 ** 6, flux_extra_int[i_theta], linestyle='-', color='m', linewidth=2)
    else :
        plt.plot(1/(100.*bande_sample)*10**6,flux_extra_int[i_theta],linestyle='-',color='r',linewidth=2)
    '''


########################################################################################################################

if Charnay == True :
    np.save('I_Charnay.npy', I_Charnay)

if D1_mean == True :

    data = np.zeros((sh_d[0],sh_d[2]))
    for i_l in range(sh_d[2]) :
        data[0,i_l] = (np.mean(data_convert[0,0,i_l,:,15]) + np.mean(data_convert[0,0,i_l,:,47]))/2.
        data[1,i_l] = (np.mean(data_convert[1,0,i_l,:,15]) + np.mean(data_convert[1,0,i_l,:,47]))/2.
        data[2,i_l] = (np.mean(data_convert[2,0,i_l,:,15]) + np.mean(data_convert[2,0,i_l,:,47]))/2.
        data[3,i_l] = (np.mean(data_convert[3,0,i_l,:,15]) + np.mean(data_convert[3,0,i_l,:,47]))/2.

    res, c_grid, i_grid = interp2olation_multi(data[0],data[1],P_comp,T_comp,x_species)
    data[6:sh_d[0]-1] = res[2:]
    for i in range(sh_d[2]) :
        data[4,i] = (1. - np.nansum(data[6:sh_d[0]-1,i]))/(ratio_HeH2 + 1.)
        data[5,i] = data[4,i]*ratio_HeH2
    h = h0

    data[sh_d[0]-1] = np.dot(M_species,data[4:sh_d[0]-1,:])

    print data, Rp, h, r_step

    I_dis, I_int = generator_1D(data,n_species,m_species,c_species,dim_bande,Rp,h,r_step,g0,cross,ind_active,k_cont,\
                 Q_cloud,P_sample,T_sample,Q_sample,bande_sample,bande_cloud,r_eff,r_cloud,rho_p,name_file,5,phi_rot,\
                 domain,ratio,directory,Tracer=False,Continuum=True,Scattering=True,Clouds=Cl,Kcorr=False,\
                 Optimal=False,Discret=True,Integral=True,Gravity=False,Molecular=True,Pressure=False,Script=False,Middle=True)

    #plt.imshow(I_dis, aspect='auto')
    #plt.colorbar()
    #plt.show()
    I_m = np.zeros((dim_bande,sh_d[2]))
    I_m[:,0] = I_dis[:,0]
    I_m[:,sh_d[2]-1] = I_dis[:,-1]
    I_m[:,1:sh_d[2]-1] = I_dis
    alt = np.zeros(sh_d[2])
    alt[0] = 0
    alt[-1] = h
    alt[1:-1] = np.linspace(r_step/2.,h-r_step/2.,sh_d[2]-2)
    for i_b in range(dim_bande) :
        alpha = 2*integrate.trapz((Rp+alt[:])*(1-I_m[i_b,:]),alt[:])
        flux_mean[i_b] = (Rp**2+alpha)/Rs**2

if Charnay == True :
    for i_b in range(dim_bande) :
        flux_extra_dis[0,i_b] = np.mean(flux_extra_dis[:,i_b])
    flux_extra_dis = flux_extra_dis[0]
    for i_b in range(dim_bande) :
        flux_extra_int[0,i_b] = np.mean(flux_extra_int[:,i_b])
    flux_extra_int = flux_extra_int[0]

    for i_b in range(dim_bande) :
        flux_dis[0,i_b] = np.mean(flux_dis[:,i_b])
    flux_dis = flux_dis[0]
    for i_b in range(dim_bande) :
        flux_int[0,i_b] = np.mean(flux_int[:,i_b])
    flux_int = flux_int[0]

########################################################################################################################

I = np.load('%sI_WASP/%s.npy'%(path,name_I))
#I = np.load('%sI/%s.npy'%(directory,name_I))

if E_O == True :

    for i in range(2):
        I_2 = np.zeros((3000,100,1))
        if i == 0 :
            I_2[:,:,0] = I[:,:,24]
        else :
            I_2[:,:,0] = I[:,:,72]

        R_eff_bar,R_eff,ratio_bar,ratR_bar,bande_bar,flux_bar,flux_n = atmospectre(I_2,bande_sample,Rs,Rp,r_step,0,\
                                                                                    False,False,True)

        if i == 0 :
            plt.plot(1 / (100. * bande_sample) * 10 ** 6, flux_n, linestyle='-', color='c', linewidth=2)
        else :
            plt.plot(1/(100.*bande_sample)*10**6,flux_n,linestyle='-',color='b',linewidth = 2)
        plt.axis([0.4, 60., 0.0189, 0.023])
        plt.show()

else :

    R_eff_bar, R_eff, ratio_bar, ratR_bar, bande_bar, flux_bar, flux_n = atmospectre(I, bande_sample, Rs, Rp, r_step,
                                                                                     0, \
                                                                                     False, False, True)


if No_extra == True :
    #pmin = np.amin(variables['p'][5,:,i_lat,i_long])
    print pmin
    da = data_convert[:,0,:,:,:]
    wh = np.where(da[0] > pmin)
    print da[0]
    print wh
    i_max = np.amax(wh[0])
    R_eff_bar,R_eff,ratio_bar,ratR_bar,bande_bar,flux_bar,flux = atmospectre(I[:,0:i_max-3,:],bande_sample,Rs,Rp,r_step,0,\
                                                                             False,False,True)
plt.figure(1)
plt.semilogx()
plt.grid(True)

if D1 == True :
    plt.plot(1/(100.*bande_sample)*10**6,flux_dis,'g',linewidth = 2,label='1D+1 discreet')
    plt.plot(1/(100.*bande_sample)*10**6,flux_int,'c',linewidth = 2,label='1D+1 integral')
if No_extra == True :
    plt.plot(1/(100.*bande_sample)*10**6,flux_n,'b',linewidth = 2,label='3D discreet extrapolate')
    plt.plot(1/(100.*bande_sample)*10**6,flux,'k',linewidth = 2,label='3D discreet')
else :
    plt.plot(1/(100.*bande_sample)*10**6,flux_n,'k',linewidth = 2,label='3D discreet')

if Charnay == True :
    plt.plot(1/(100.*bande_sample)*10**6,flux_extra_dis,'r',linewidth = 2,label='3D+1 discreet')
    #plt.plot(1/(100.*bande_sample)*10**6,flux_extra_int,'r',linewidth = 2,label='3D+1 integral')
if D1_mean == True :
    plt.plot(1/(100.*bande_sample)*10**6,flux_mean,'g',linewidth = 2,label='1D mean profile')
plt.ylabel('Relative flux')
plt.xlabel('Wavelenght (micron)')
plt.draw()

plt.figure(2)
plt.semilogx()
plt.grid(True)

#plt.plot(1/(100.*bande_sample)*10**6,(flux_n - flux_dis)*1.e+6,'g',linewidth = 2,label='1D+1 discreet')
#plt.plot(1/(100.*bande_sample)*10**6,(flux_n - flux_int)*1.e+6,'c',linewidth = 2,label='1D+1 integral')
#plt.plot(1/(100.*bande_sample)*10**6,(flux - flux_extra_int)*1.e+6,'b',linewidth = 2,label='3D discreet extrapolate')
if D1_mean == True :
    plt.plot(1/(100.*bande_sample)*10**6,np.abs(flux_n - flux_mean)*1.e+6,'k',linewidth = 2,label='3D discreet')
if Charnay == True :
    plt.plot(1/(100.*bande_sample)*10**6,np.abs(flux_n - flux_extra_dis)*1.e+6,'k',linewidth = 2,label='3D discreet')
if No_extra == True :
    plt.plot(1/(100.*bande_sample)*10**6,np.abs(flux_n - flux)*1.e+6,'k',linewidth = 2,label='Extrapolation')
#plt.plot(1/(100.*bande_sample)*10**6,(flux_n - flux_extra_dis)*1.e+6,'y',linewidth = 2,label='3D+1 discreet')
#plt.plot(1/(100.*bande_sample)*10**6,(flux_n - flux_extra_int)*1.e+6,'r',linewidth = 2,label='3D+1 integral')
if Noise == True :
    detection = JWST()
    class star :
        def __init__(self):
            self.radius = Rs
            self.temperature = Ts
            self.distance = d_al
    bande_sample = np.delete(bande_sample,[0])
    int_lambda = np.zeros((2,bande_sample.size))
    bande_sample = np.sort(bande_sample)
    resolution = 'low'
    if resolution == '' :
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
        print int_lambda
    else :
        int_lambda = np.sort(10000./bande_sample[::-1])

    noise = stellar_noise(star(),detection,int_lambda,resolution)
    noise = noise[::-1]
    plt.plot(1/(100.*bande_sample)*10**6,noise*10**6,'m',linewidth=2,label='JWST')

plt.ylabel('Relative error flux (ppm)')
plt.xlabel('Wavelenght (micron)')
plt.show()