import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
sys.path.append('/Users/caldas/Pytmosph3R/PyRouts/')
from pyfunction import *
from pyconstant import *

#directory = '/Users/caldas/Desktop/Pytmosph3R/Tools/'
#name_data = 'GJ1214b_data_convert_120006448.npy'
directory = '/Users/caldas/'
name_data = 'Desktop/Heterogeneities/GJ1214b/Para_GJ1214b/Parameters/GJ1214b_data_convert_6961x64x48_lo0.00_0.5mu_0.5mu.npy'
data_convert = np.load('%s/%s'%(directory,name_data))
data_convert = data_convert[:,:,:82,:,:]
sh = np.shape(data_convert)
parameters = np.array(['P','T','gen_cond1','gen_cond2','H2','He','H2O','CH4','N2','NH3','CO','CO2','M'])
#stud = np.array(['P','T','gen_cond1','gen_cond2','CH4','CO2','CO','H2O','NH3','M'])
stud = np.array(['P','T','CH4','H2O','NH3','CO2','CO'])

mi = np.array([0,300,-6.8,-1.7,-9,-10,-10])

i_lat = 24
factor = 10.

n_layers = 80
r_step = 69615.
#r_step = 120000.
h = n_layers*r_step
param = 20000.
Rp = 1.7e+7
#Rp = 15.*R_T
Rp_m = Rp/factor
dim = int((Rp_m+1.10*h)/param)

reso_long = 64
long_step = 2*np.pi/np.float(reso_long)

data_map = data_convert[:,0,:,i_lat,:]
print np.shape(data_map)
p_min = np.amin(data_map[0,:,24])
p_min = 3.e-1

map = np.zeros((stud.size,2*dim+1,2*dim+1))
map_ref = np.zeros((2*dim+1,2*dim+1))
map2 = np.zeros((2*dim+1,2*dim+1))
i_X = 0

bar = ProgressBar(2*dim+1,'Progression : ')

for i_x in range(-dim,dim+1) :
    i_Y = 0
    for i_y in range(-dim,dim+1) :
        x = i_x*param
        y = i_y*param

        if np.sqrt((x**2 + y**2)) >= Rp_m and np.sqrt((x**2 + y**2)) <= Rp_m + h :

            i_z = np.int((np.sqrt((x**2 + y**2))-Rp_m)/r_step)+1
            if np.sqrt((x**2 + y**2))-Rp_m < r_step/2. :
                i_z = 0
            if np.sqrt((x**2 + y**2))-Rp_m > h-r_step/2. :
                i_z = n_layers+1
            long = np.abs(np.arccos(np.abs(x)/np.sqrt((x**2 + y**2))))

            if long%long_step < long_step/2. :
                modif = 0
            else :
                modif = 1

            if x < 0 :
                i_long = np.int(reso_long/2 - (i_y+0.1)/np.abs(i_y+0.1)*np.int(long/long_step + modif))
            else :
                if y > 0 :
                    i_long = (np.int(long/long_step) + modif)
                else :
                    i_long = reso_long - np.int(long/long_step) - modif

            for i_st in range(stud.size) :

                i_stud, = np.where(parameters == stud[i_st])
                i_stud = np.int(i_stud[0])
                if stud[i_st] == 'T' or stud[i_st] == 'M' :
                    map[i_st,i_Y,i_X] = data_map[i_stud,i_z,i_long]
                else :
                    if stud[i_st] == 'P' :
                        if data_map[0,i_z,i_long] >= p_min :
                            map[i_st,i_Y,i_X] = np.log10(data_map[i_stud,i_z,i_long])
                        else :
                            map[i_st,i_Y,i_X] = 'nan'
                        map_ref[i_Y,i_X] = np.log10(data_map[0,i_z,i_long])
                    else :
                        if stud[i_st] == 'gen_cond1' or stud[i_st] == 'gen_cond2' :
                            if data_map[i_stud,i_z,i_long] <= 1.e-7 :
                                map[i_st, i_Y, i_X] = 'nan'
                            else :
                                map[i_st,i_Y,i_X] = np.log10(data_map[i_stud,i_z,i_long])
                        else :
                            if data_map[i_stud,i_z,i_long] <= 1.e-10 :
                                map[i_st, i_Y, i_X] = -10
                            else :
                                map[i_st, i_Y, i_X] = np.log10(data_map[i_stud, i_z, i_long])
                #if i_st == 0 :
                #    map2[i_Y,i_X] = i_long
        else :

            for i_st in range(stud.size) :

                map[i_st,i_Y,i_X] = 'nan'

        i_Y += 1
    i_X += 1
    bar.animate(i_X)

for i_st in range(stud.size) :

    print 'hola'
    x = np.linspace(-dim*param,dim*param,2*dim+1)
    y = x
    X,Y = np.meshgrid(x,y)
    lev = np.array([np.log10(p_min)])
    x_cir,y_cir = create_circle(x,Rp_m)
    x_cir_r,y_cir_r = create_circle(x,Rp_m+h)

    plt.figure(figsize=(12, 9))
    plt.grid(True)
    plt.imshow(map[i_st],aspect='auto',origin='lower',extent=np.array([-dim*param/1000.,dim*param/1000.,-dim*param/1000.,dim*param/1000.]),cmap=cm.jet,vmin=mi[i_st],vmax = 650)
    plt.colorbar()
    print 'hola'
    '''
    if stud[i_st] != 'P' :
        CS = plt.contour(X/1000.,Y/1000.,map_ref,levels=lev,colors='w',linewidths=2)
        #plt.clabel(CS,frontsize = 3, inline = 0)
    else :
        CS = plt.contour(X/1000.,Y/1000.,map_ref,levels=lev,colors='w',linewidths=0)
    '''
    plt.plot(y_cir/1000.,x_cir/1000.,'-k',linewidth=2)
    plt.plot(y_cir_r/1000.,x_cir_r/1000.,linestyle='--',color='k',linewidth=3)
    I_min = np.load('I_min.npy')
    I_max = np.load('I_max.npy')
    CS = plt.contour(X / 1000., Y / 1000., I_min, levels=np.array([0.05]), colors='k', linewidths=0)
    CS = plt.contour(X / 1000., Y / 1000., I_max, levels=np.array([0.95]), colors='k', linewidths=0)
    plt.show()

    #plt.grid(True)
    #plt.imshow(map2,aspect='auto',origin='lower',interpolation='none',extent=np.array([-dim*step,dim*step,-dim*step,dim*step]))
    #plt.colorbar()
    #plt.show()

