import numpy as np
import matplotlib as plt
import pickle as pkl

#data = pkl.load(open('nest_out.pickle','rb'))
data = pkl.load(open('GJ1214b_Est/nest_out.pickle','rb'))

data.keys()

sh = np.shape(data['solutions'])

y = np.array(['radius','T','log_H2O','log_NH3','log_CO2','log_CO','log_N2','log_CH4'])
#y = np.array(['radius','T','log_H2O','log_NH3','log_CO2','log_CO','log_CH4'])

for i in range(8) :
    array = data['solutions'][0]['fit_params'][y[i]]['value']
    error_plus = data['solutions'][0]['fit_params'][y[i]]['sigma_p']
    error_moins = data['solutions'][0]['fit_params'][y[i]]['sigma_m']

    print y[i], array, error_plus, error_moins