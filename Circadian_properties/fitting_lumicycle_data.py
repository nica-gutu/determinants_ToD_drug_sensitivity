#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 11:19:44 2024

@author: nicagutu
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

#function to fit
def expsin(time, A0, lam, T, phase):
    return A0*np.exp(-lam*time)*np.sin(2*np.pi/T*time+phase)

path = 'output/'

path2 = 'Raw_Data/'
file = 'lumicycle_data_breast_U2OS_5to142h.xlsx'
sheet_name = 'detrended'

df = pd.read_excel(path2+file, sheet_name=sheet_name)
df = df.fillna(0)

time = df['Time']
time_fit = np.linspace(0, max(df['Time']), len(df['Time']))

initial_guess = [1.0, 0.1, 1.0, 0.0]  # [A0, lam, T, phase]
bounds = ([0, 0, 0, -2*np.pi], [np.inf, np.inf, np.inf, 2*np.pi])  # Bounds for all parameters

dt = 1/6 
lowT = 16
highT = 50
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#to save parameters from fitting
Ampl_param = []
Period_param = []
decay_param = []
#to save fitting errors
Ampl_error = []
Period_error = []
decay_error = []

label = []

count = 0
for col in df.columns[1:len(df.columns)]:
    signal = df[col]
    trend = wAn.sinc_smooth(signal, T_c=10)
    time_fit = np.linspace(0, max(time), len(time))

    try:
        popt, pcov = curve_fit(expsin, time_fit, trend, p0=initial_guess, bounds=bounds, method='trf', maxfev=100000)
        param_errors = np.sqrt(np.diag(pcov))
 
        # Calculating R-squared
        y_fit = expsin(time_fit, *popt)
        y_obs = signal
        y_mean = np.mean(y_obs)
        SST = np.sum((y_obs-y_mean)**2)
        SSE = np.sum((y_obs-y_fit)**2)
        R_squared = 1-(SSE/SST)

        if R_squared <= 0.5: #excluding bad fitting with R2 smaller than 0.5
            count += 1            
            print('R-squared: ', R_squared)
            print('Parameters: ', popt)
        
            plt.figure(figsize=(10, 6))
            plt.plot(time_fit, signal, 'ko', label='Original Data')
            plt.plot(time_fit, trend, '-', label='Trend')
            plt.plot(time_fit, expsin(time_fit, *popt), 'r-', label='Fitted Function')
            plt.xlabel('Time')
            plt.ylabel('df[col]')
            plt.title(str(col))
            plt.legend()
            plt.show()
        else:
            Ampl_param.append(popt[0])
            Period_param.append(popt[2])
            decay_param.append(popt[1])
            Ampl_error.append(param_errors[0])
            Period_error.append(param_errors[2])
            decay_error.append(param_errors[1])
            label.append(col)


    except RuntimeError as e:
        print(f"Error occurred for column {col}: {e}")
        continue

# print(Ampl_param)
# print(Period_param)
# print(decay_param)
# print(label)

print(count)
    
df_properties = pd.DataFrame(columns=['Label', 'Amplitude', 'Error_A', 'Period', 'Error_T', 'decay', 'Error_d'])
df_properties['Label'] = label
df_properties['Amplitude'] = Ampl_param
df_properties['Error_A'] = Ampl_error
df_properties['Period'] = Period_param
df_properties['Error_T'] = Period_error
df_properties['decay'] = decay_param
df_properties['Error_d'] = decay_error
# print(df_properties)

df_properties.to_excel(path2+'Fitting_lumicycle_data.xlsx', index=None)

#%%Group by labels
# cleaned_labels = []
# for col in label:
#     cleaned_label = col.split('_', 1)[1]
#     cleaned_label = cleaned_label.rsplit('_', 1)[0]
#     cleaned_labels.append(cleaned_label)
    
# df_properties.columns = cleaned_labels    

# averaged_data = df_properties.groupby(df_properties.columns, axis=1).mean()
# print(averaged_data)

# (averaged_data.T).to_excel(path+'Circadian_properties_new_fitting.xlsx')



