#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:10:32 2024

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

def GR_2(c, Sm, SC50, h):
    return 2**(1-(Sm*c**h)/(SC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

path ='Raw_Data/'
path2 = 'output/'
file = '041922_Dose_Response'

data = pd.read_excel(path + file + '.xlsx', header=0, sheet_name='CellCount')
indexs = data.index.values
eval_time = [48, 72, 96, 120]#, len(data)-1]

cols = ['B','C','D','E','F','G'] 

list_drugs = ['5-FU','Cisplatin','Doxorubicin']

concentrations_5FU = [100, 50, 25, 12.5, 6.25, 3.13, 1.56, 0.78, 0.39]  # 5-FU
concentrations_Cisp = [70, 35, 17.5, 8.755, 4.38, 2.19, 1.09, 0.55, 0.27]  # Cisplatin
concentrations_Dox = [1, 0.5, 0.25, 0.13, 0.06, 0.03, 0.02, 0.01, 3.91e-3]  # Doxorubicin
concentrations_all = (concentrations_5FU, concentrations_Cisp, concentrations_Dox)

inf_coeffs = np.zeros((len(eval_time), len(list_drugs)))
EC50_coeffs = np.zeros((len(eval_time), len(list_drugs)))
h_coeffs = np.zeros((len(eval_time), len(list_drugs)))

#%%Survival curves
for a, b, c in zip([0,2,4], list_drugs, range(3)):
    drug = b

    plt.figure(figsize=(12,10))

    for ii in range(len(eval_time)):
        tp = eval_time[ii]
        
        concentrations = concentrations_all[c]
        concs2 = np.arange(concentrations[-1], concentrations[0], step=0.001)
        
        i = cols[a]
        j = cols[a+1]
        pattern1 = re.compile(rf'(?:{i})[2-9]|(?:{i})1[01]')
        pattern2 = re.compile(rf'(?:{j})[2-9]|(?:{j})1[01]')  
    
        selected_columns1 = data.filter(regex=pattern1)
        selected_columns2 = data.filter(regex=pattern2)
        
        control1 = selected_columns1.pop(selected_columns1.columns[0])
        control2 = selected_columns2.pop(selected_columns2.columns[0])
        mean_control = (control1+control2)/2
        control_std_dev = np.std([control1, control2], axis=0)
            
        DRC = []
        error_DRC = []
        for n,m in zip(selected_columns1, selected_columns2):
            mean = (selected_columns1[n]+selected_columns2[m])/2
            std_dev = np.std([selected_columns1[n], selected_columns2[m]], axis=0)
            
            point = (mean[tp]/mean_control[tp])
            error = np.abs(point)*(np.sqrt((std_dev[tp]/mean[tp])**2+(control_std_dev[tp]/mean_control[tp])**2))       
            
            DRC.append(point)
            error_DRC.append(error)
        
        DRC = DRC/DRC[-1]
        
        # print(DRC)
        # print(concentrations)
        params, cov = curve_fit(GR, concentrations, DRC, p0=[0.2, 0.5, 2])
        inf_coeffs[ii,c] = params[0]
        EC50_coeffs[ii,c] = params[1]
        h_coeffs[ii,c] = params[2]

        plt.errorbar(concentrations, DRC, yerr=error_DRC, fmt='o', markersize=12, label=str(eval_time[ii]))
        plt.plot(concs2, GR(concs2, *params)) ###or GR
    plt.xscale('log')
    plt.xlabel('Dose')
    plt.ylabel('Survival fraction')
    plt.legend(loc='best', title='Evaluation')
    plt.title(str(b))
    # plt.savefig(path2+'Survival_time_eval_'+str(b)+'.svg')
    plt.show()

df_inf = pd.DataFrame(inf_coeffs, columns=list_drugs, index=eval_time)                
df_EC50 = pd.DataFrame(EC50_coeffs, columns=list_drugs, index=eval_time)                
df_h = pd.DataFrame(h_coeffs, columns=list_drugs, index=eval_time)      

#%%Simulated ToD response       

#Simulation parameters
toD = np.arange(0, 24.1, step=.1) #time of day
tf = 120
t = np.linspace(0, tf, tf)

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5
colors2 = plt.cm.tab20b(range(len(toD)))
 
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate

#Drug parameters
Sm = 0 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect

step = 24
toE = np.arange(eval_time[0], eval_time[-1]+step, step=step)

list_drugs = ['Cisplatin']
for a in list_drugs:
    plt.figure(figsize = (12,10))
    EC50 = df_EC50[a][120]
    kl = df_inf[a][120]

    for b in eval_time:
        ratio = []
        h = df_h[a][b]    
        c2 = EC50*circmod_c(modul_strgth, Amp, toD, T, phi)

        for j in range(len(c2)):
            control = growthcurve(1, b, k, kl, LC50, Sm, 0, SC50, h)
            ratio.append(growthcurve(1, b, k, kl, LC50, Sm, c2[j], SC50, h)/control)        

        plt.plot(toD, ratio/ratio[0], linewidth=4, label=str(b)+'h')
    plt.ylabel('ToD response')
    plt.xticks([0,6,12,18,24])
    plt.xlabel('Time of day [hours]')
    plt.title(str(a))
    plt.legend(loc='best', title='Evaluation time')
    # plt.savefig(path2+'ToD_expected_from_exp_'+str(a)+'_eval_time.svg')
    plt.show()

















      
                
