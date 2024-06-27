#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:50:18 2024

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
from scipy.optimize import least_squares

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

def exp(x, a, b, c):
    return a*np.exp(b*np.array(x))+c

def residual(params, x, y):
    return y-exp(x, *params)

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
file = '081122_Drug_Stability'

conditions2 = ['B', 'C', 'D'] 
conditions3 = ['E', 'F', 'G'] 
list_drugs = ['5FU', 'Dox', 'Cisplatin']
list_cond = ['48h', '24h', '0h']

df_EC50 = pd.DataFrame(columns=list_drugs, index=list_cond)
df_EC50_errors = pd.DataFrame(columns=list_drugs, index=list_cond)

df_Einf = pd.DataFrame(columns=list_drugs, index=list_cond)
df_h = pd.DataFrame(columns=list_drugs, index=list_cond)

concentrations_5FU = [100, 50, 25, 12.5, 6.25, 3.13, 1.56, 0.78, 0.39]  # 5-FU
concentrations_Cisp = [70, 35, 17.5, 8.755, 4.38, 2.19, 1.09, 0.55, 0.27]  # Cisplatin
concentrations_Dox = [1, 0.5, 0.25, 0.13, 0.06, 0.03, 0.02, 0.01, 3.91e-3]  # Doxorubicin
concentrations_all = {'5FU':concentrations_5FU, 'Dox':concentrations_Dox, 'Cisplatin': concentrations_Cisp}

# concentrations = [1, 0.5, 0.25, 0.13, 0.06, 0.03, 0.02, 0.01, 3.91e-3]
# concs2 = np.arange(1e-3, 1, step=1e-3)

#%%Survival curves
for a, b, c in zip([1,2,3], [2,3,1], list_drugs):
    data1 = pd.read_excel(path + file + '.xlsx', header=0, sheet_name=f'Plate{a}_CellCount')
    data2 = pd.read_excel(path + file + '.xlsx', header=0, sheet_name=f'Plate{b}_CellCount')
    
    time = data1.index.values
    last_tp = len(data2)-1

    plt.figure(figsize=(10,8))
    EC50_coeffs = []
    EC50_errors = []
    Einf_coeffs = []
    h_coeffs = []
    
    concentrations = (concentrations_all[c])
    concs2 = np.arange(concentrations[-1], concentrations[0], step=0.0001)
    
    for i, j, k in zip(conditions2, conditions3, list_cond):
        pattern1 = re.compile(rf'(?:{i})[2-9]|(?:{i})1[01]')
        pattern2 = re.compile(rf'(?:{j})[2-9]|(?:{j})1[01]')  

        selected_columns1 = data1.filter(regex=pattern1)
        selected_columns2 = data2.filter(regex=pattern2)
                                
        control1 = selected_columns1.pop(selected_columns1.columns[0])
        control2 = selected_columns2.pop(selected_columns2.columns[0])
        mean_control = (control1+control2)/2
        control_std_dev = np.std([control1, control2], axis=0)
        
        DRC = []
        error_DRC = []
        for n,m in zip(selected_columns1, selected_columns2):
            mean = (selected_columns1[n]+selected_columns2[m])/2
            std_dev = np.std([selected_columns1[n], selected_columns2[m]], axis=0)
            point = (mean[last_tp]/mean_control[last_tp])
            error = np.abs(point)*(np.sqrt((std_dev[last_tp]/mean[last_tp])**2+(control_std_dev[last_tp]/mean_control[last_tp])**2))       
            DRC.append(point)
            error_DRC.append(error)
        DRC = DRC/DRC[-1]
        plt.errorbar(concentrations, DRC, yerr=error_DRC, fmt='o', label=str(k))
        
        # params, cov = curve_fit(GR, concentrations, DRC, p0=[0.2, 0.5, 2])
        params, cov = curve_fit(GR_2, concentrations, DRC, p0=[0.2, 0.5, 2], maxfev=20000, method='lm')
        errors = np.sqrt(np.diag(cov))
        
        Einf_coeffs.append(params[0])
        EC50_coeffs.append(params[1])
        h_coeffs.append(params[2])
        EC50_errors.append(errors[1])
    
        plt.plot(concs2, GR_2(concs2, *params)) ###or GR
    
    plt.ylabel('Survival')
    plt.xlabel('Dose $\mu$M')        
    plt.xscale('log')
    plt.legend(loc='best')
    plt.title(str(c))
    # plt.savefig(path2+'Survival_experimental_'+str(c)+'_drug_age.svg')
    plt.show()
    
    df_EC50[c] = EC50_coeffs
    df_EC50_errors[c] = EC50_errors

    df_Einf[c] = Einf_coeffs
    df_h[c] = h_coeffs


#%%EC50 

list_cond2 = [0,24,48]

decay_rates = []

plt.figure(figsize=(10,8))
for i in df_EC50:
    delta_A = df_EC50_errors[i]
    delta_B = df_EC50_errors[i][2]
    ratio = df_EC50[i]/df_EC50[i][2] 
    propagated_error = np.abs(ratio)*np.sqrt((delta_A/df_EC50[i])**2+(delta_B/df_EC50[i][2])**2)
    
    if i !='5FU':
        plt.errorbar(list_cond2, ratio, yerr=propagated_error, fmt='-o', markersize=12, label=str(i))
    if i == '5FU':
        plt.plot(list_cond2, ratio, '-o', markersize=12, label=str(i))
            
    initial_guess = (1, -0.01, 1)
    result = least_squares(residual, initial_guess, args=(list_cond2, ratio))
    popt = result.x
    print(popt)
    decay_rates.append(abs(popt[1]))
    
    x_fit = np.linspace(min(list_cond2), max(list_cond2), 100)
    plt.plot(x_fit, exp(x_fit, *popt), color='black', linestyle='--', label='Fitting '+str(i))
    
plt.ylabel('Relative EC50 to 0h')
plt.xlabel('Drug age [hours]')
plt.legend(loc='best')
# plt.ylim([0, 5])
plt.gca().invert_xaxis()
plt.xticks([0,24,48],[48,24,0])
# plt.savefig(path2+'Rel_EC50_drug_age.svg')
plt.show()

#%%Simulated ToD response
      
#Simulation parameters
toD = np.arange(0, 25, step=1) #time of day
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

plt.figure(figsize = (12,10))
for ii in range(1,3):
    a = list_drugs[ii]
    print(a, decay_rates[ii])
    EC50 = 1#df_EC50[a][-1]
    c2 = EC50*circmod_c(modul_strgth, Amp, toD, T, phi)*np.exp(-(decay_rates[ii])*toD)
    Sm = df_Einf[a][-1]
    h = df_h[a][-1]
    ratio = []
    for j in range(len(c2)):
        control = growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1]
        treated = growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]
        ratio.append(treated/control)
    plt.plot(toD, ratio/ratio[0], linewidth=4, label=str(a)+': decay {:.3f}'.format(decay_rates[ii])+'h')
c2 = circmod_c(modul_strgth, Amp, toD, T, phi)#*np.exp(-(decay_mean[ii])*toD)
Sm = 1#df_Einf[a][-1]
h = 1#df_h[a][-1]
ratio = []
for j in range(len(c2)):
    control = growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1]
    treated = growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]
    ratio.append(treated/control)
plt.plot(toD, ratio/ratio[0], linewidth=4, label=r'Decay at $\infty$ h')
plt.xticks([0,6,12,18,24])
plt.ylabel('ToD response')
plt.xlabel('Time of day [hours]')
plt.legend(loc='best')
# plt.savefig(path2+'ToD_expected_from_exp_Dox_Cispl_norm_drug_age.svg')
plt.show()








