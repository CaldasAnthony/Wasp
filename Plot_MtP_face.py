import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
sys.path.append('/Users/caldas/Pytmosph3R/ComnPyRouts/')
from pyfunction import *
from pyconstant import *

directory = '/Users/caldas/'
name_data = 'Desktop/Heterogeneities/GJ1214b/Para_GJ1214b/Parameters/GJ1214b_data_convert_6961x64x48_lo0.00_0.5mu_0.5mu.npy'
data_convert = np.load('%s%s'%(directory,name_data))
data_convert = data_convert[:,:,:82,:,:]
sh = np.shape(data_convert)
parameters = np.array(['P','T','gen_cond1','gen_cond2','H2','He','H2O','CH4','N2','NH3','CO','CO2','M'])
stud = np.array(['P','T','gen_cond1','gen_cond2','CH4','CO2','CO','H2O','NH3','M'])

i_right = 16
i_left = 48
factor = 10.
theta_number = 96

#n_layers = 100
#r_step = 52120.
n_layers = 80
r_step = 69610.
h = n_layers*r_step
param = 20000.
Rp = 1.7e+7
Rs = 0.206470165349*R_S
Rp_m = Rp/factor
dim = int((Rp_m+1.10*h)/param)
print dim

reso_lat = 48
lat_step = np.pi/np.float(reso_lat)

data_map = np.zeros((sh[0],2,sh[2],sh[3]))
data_map[:,0] = data_convert[:,0,:,:,i_right]
data_map[:,1] = data_convert[:,0,:,:,i_left]
p_min = np.amin(data_map[0,0,:,24])
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
            lat = np.abs(np.arccos(np.abs(x)/np.sqrt((x**2 + y**2))))

            if lat%lat_step < lat_step/2. :
                if i_y >=0 :
                    modif = 0
                else :
                    modif = 0
            else :
                if i_y >= 0 :
                    modif = 1
                else :
                    modif = -1
            i_lat = np.int(reso_lat/2 + (i_y+0.1)/np.abs(i_y+0.1)*np.int(lat/lat_step) + modif)

            if i_x <= 0 :
                i_long = 1
            else :
                i_long = 0

            for i_st in range(stud.size) :

                i_stud, = np.where(parameters == stud[i_st])
                i_stud = np.int(i_stud[0])
                if stud[i_st] == 'T' or stud[i_st] == 'M' :
                    map[i_st,i_Y,i_X] = data_map[i_stud,i_long,i_z,i_lat]
                else :
                    if stud[i_st] == 'P' :
                        if data_map[0,i_long,i_z,i_lat] >= p_min :
                            map[i_st,i_Y,i_X] = np.log10(data_map[i_stud,i_long,i_z,i_lat])
                        else :
                            #map[i_st,i_Y,i_X] = 'nan'
                            map[i_st, i_Y, i_X] = np.log10(data_map[i_stud, i_long, i_z, i_lat])
                        map_ref[i_Y,i_X] = np.log10(data_map[0,i_long,i_z,i_lat])
                    else :
                        if data_map[i_stud,i_long,i_z,i_lat] <= 1.e-6 :
                            map[i_st, i_Y, i_X] = 'nan'
                        else :
                            map[i_st,i_Y,i_X] = np.log10(data_map[i_stud,i_long,i_z,i_lat])

                #if i_st == 0 :
                #    map2[i_Y,i_X] = i_lat
        else :

            for i_st in range(stud.size) :

                map[i_st,i_Y,i_X] = 'nan'

        i_Y += 1
    i_X += 1
    bar.animate(i_X)

lev = np.array([np.log10(p_min)])

for i_st in range(stud.size) :

    x_cir,y_cir = create_circle(x,Rp_m)
    x_cir_r,y_cir_r = create_circle(x,Rp_m+h)

    plt.figure(figsize=(12, 9))
    plt.grid(True)
    plt.imshow(map[i_st],aspect='auto',origin='lower',extent=np.array([-dim*param/1000.,dim*param/1000.,-dim*param/1000.,dim*param/1000.]),cmap=cm.jet)
    plt.colorbar()
    if stud[i_st] != 'P' :
        CS = plt.contour(X/1000.,Y/1000.,map_ref,levels=lev,colors='w',linewidths=2)
        #plt.clabel(CS,frontsize = 3, inline = 0)
    plt.plot(y_cir/1000.,x_cir/1000.,'-k',linewidth=2)
    plt.plot(y_cir_r/1000.,x_cir_r/1000.,linestyle='--',color='k',linewidth=3)
    plt.show()

    #plt.grid(True)
    #plt.imshow(map2,aspect='auto',origin='lower',interpolation='none',extent=np.array([-dim*step,dim*step,-dim*step,dim*step]))
    #plt.colorbar()
    #plt.show()
