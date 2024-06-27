#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 12:08:45 2024

@author: nicagutu
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import numpy as np

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

def linear_func(x, m, b):
    return m*x+b

def exponential_func(x, a, b):
    return a*np.exp(b*x)

path ='Raw_Data/'
path2 = 'output/'
file = 'Circadian_properties_lumicycle'

data = pd.read_excel(path + file + '.xlsx', header=0, sheet_name='Maximum_ToD')

high_corr_ampl = ['Alpelisib']#, '5-FU', 'Cisplatin', 'Torin2']
high_corr_per = ['Alpelisib'] #,'Alisertib', 'Doxorubicin',
high_corr_decay = ['Alpelisib']

#heatmap with all correlation coefficients
correlation_matrix = data.corr()
plt.figure(figsize=(14, 12))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", annot_kws={"size": 15},
            cbar_kws={"label": "Correlation Coefficient"})
plt.show()

celllines = data['Cellline']
print(celllines)
periodBmal1_2d = data[data.columns[1]]
periodPer_2d = data[data.columns[2]]
periodBmal1_5d = data[data.columns[3]]
periodPer_5d = data[data.columns[4]]
amplBmal1 = data[data.columns[5]]
amplPer = data[data.columns[6]]
decayBmal1 = data[data.columns[7]]
decayPer = data[data.columns[8]]

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

#%%Amplitude
plt.figure(figsize=(12,10)) 

for col in high_corr_ampl:
    max_tod = data[col].dropna()
    max_tod = (max_tod-min(max_tod))/(max(max_tod)-min(max_tod))
    amplPer_filtered = amplPer[max_tod.index]
        
    for j in max_tod.index:
        plt.plot(amplPer_filtered[j], max_tod[j], 'o', markersize=10, color=colors[j], label=celllines[j])
    
    popt, pcov = curve_fit(linear_func, amplPer_filtered, max_tod)
    corr_coef, p_value = pearsonr(amplPer_filtered, max_tod)

plt.xlabel('Amplitude Per [a.u.]')
plt.ylabel('Maximum relative ToD response')
# plt.legend(loc='best')
# plt.savefig(path2+'Amplitude_max_ToD_Alisertib.svg')
plt.show()


#%%Period
for col in high_corr_per:
    max_tod = data[col].dropna()
    periodBmal1_filtered = periodBmal1_5d[max_tod.index] 
        
    plt.figure(figsize=(12,10))
    for j in max_tod.index:
        plt.plot(periodBmal1_filtered[j], max_tod[j], 'o', markersize=10, label=celllines[j])

    popt, pcov = curve_fit(linear_func, periodBmal1_filtered, max_tod)
    corr_coef, p_value = pearsonr(periodBmal1_filtered, max_tod)
    plt.plot(periodBmal1_filtered, linear_func(periodBmal1_filtered, *popt), '-k', label='Corr.={:.2f}, p-value={:.3f}'.format(corr_coef, p_value))
           
    plt.xlabel('Period Bmal1 [hours]')
    plt.ylabel('Maximum relative ToD response')
    plt.legend(loc='best')
    # plt.savefig(path2+'Period_max_ToD_Alisertib.svg')
    plt.show()


#%%Decay
for col in high_corr_decay:    
    max_tod = data[col].dropna()
    decayPer_filtered = decayPer[max_tod.index]

    plt.figure(figsize=(12,10))    
    for j in max_tod.index:
        plt.plot(decayPer_filtered[j], max_tod[j], 'o', markersize=10, label=celllines[j])

    popt, pcov = curve_fit(linear_func, decayPer_filtered, max_tod)
    corr_coef, p_value = pearsonr(decayPer_filtered, max_tod)
    plt.plot(decayPer_filtered, linear_func(decayPer_filtered, *popt), '-k', label='Linear Fit (corr.={:.2f}, p-value={:.3f})'.format(corr_coef, p_value))
    
    plt.xlabel('Decay Per [1/hours]')
    plt.ylabel('Maximum relative ToD response')
    plt.legend(loc='best')
    # plt.savefig(path2+'Decay_max_ToD_Alisertib.svg')
    plt.show()

#%% 2nd fitting

data = pd.read_excel(path + file + '.xlsx', header=0, sheet_name='Maximum_ToD2')

high_corr_ampl = ['5FU_confl', 'Torin2_confl', 'Alpelisib_confl', 'Dox_confl', 'Adavosertib_confl', 'Alisertib_confl']
high_corr_per = ['Cispl_confl'] 
high_corr_decay = ['Torin2_cell', 'Dox_cell']

correlation_matrix = data.corr()
plt.figure(figsize=(16, 14))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", annot_kws={"size": 10},
            cbar_kws={"label": "Correlation Coefficient"})
plt.show()

celllines = data['Cellline']
print(celllines)

amplBmal1 = data[data.columns[1]]
amplBmal1_e = data[data.columns[2]]
periodBmal1 = data[data.columns[3]]
periodBmal1_e = data[data.columns[4]]
decayBmal1 = data[data.columns[5]]
decayBmal1_e = data[data.columns[6]]
amplPer = data[data.columns[7]]
amplPer_e = data[data.columns[8]]
periodPer = data[data.columns[9]]
periodPer_e = data[data.columns[10]]
decayPer = data[data.columns[11]]
decayPer_e = data[data.columns[12]]

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

#%%Amplitude

for col in high_corr_ampl:
    max_tod = data[col].dropna()
    amplBmal1_filtered = amplBmal1[max_tod.index].dropna()   
    max_tod = max_tod[amplBmal1_filtered.index]

    nums = (len(max_tod.index))
    
    corr_coef, p_value = pearsonr(amplBmal1_filtered, max_tod)
    
    if nums >= 10 and p_value < 0.05:
        print(col, nums, p_value)
        plt.figure(figsize=(12,10)) 
        
        for j in max_tod.index:
            plt.errorbar(amplBmal1_filtered[j], max_tod[j], xerr=amplBmal1_e[j], fmt='o', markersize=10, label=celllines[j])
        
        popt, pcov = curve_fit(linear_func, amplBmal1_filtered, max_tod)
        plt.plot(amplBmal1_filtered, linear_func(amplBmal1_filtered, *popt), '-k', label=str(col)+'\nCorr.={:.2f}, p-value={:.3f}'.format(corr_coef, p_value))
    
        plt.xlabel('Amplitude Bmal1 [a.u.]')
        plt.ylabel('Maximum relative ToD response')
        plt.legend(loc='best')
        plt.savefig(path2+'Amplitude_max_ToD_'+str(col)+'.svg')
        plt.show()


#%%Period
for col in high_corr_per:
    max_tod = data[col].dropna()
    periodPer_filtered = periodPer[max_tod.index].dropna()   
    max_tod = max_tod[periodPer_filtered.index]

    nums = (len(max_tod.index))
    
    corr_coef, p_value = pearsonr(periodPer_filtered, max_tod)
    
    if nums >= 10 and p_value < 0.05:
        print(col, nums, p_value)
        plt.figure(figsize=(12,10)) 
        
        for j in max_tod.index:
            plt.errorbar(periodPer_filtered[j], max_tod[j], xerr=periodPer_e[j], fmt='o', markersize=10, label=celllines[j])
    
        popt, pcov = curve_fit(linear_func, periodPer_filtered, max_tod)
        plt.plot(periodPer_filtered, linear_func(periodPer_filtered, *popt), '-k', label=str(col)+'\nCorr.={:.2f}, p-value={:.3f}'.format(corr_coef, p_value))
    
        plt.xlabel('Period Per2 [hours]')
        plt.ylabel('Maximum relative ToD response')
        plt.legend(loc='best')
        plt.savefig(path2+'Period_max_ToD_'+str(col)+'.svg')
        plt.show()



#%%Decay
for col in high_corr_decay:
    max_tod = data[col].dropna()
    decayBmal1_filtered = decayBmal1[max_tod.index].dropna()   
    max_tod = max_tod[decayBmal1_filtered.index]

    nums = (len(max_tod.index))
    
    corr_coef, p_value = pearsonr(decayBmal1_filtered, max_tod)
    
    if nums >= 10 and p_value < 0.05:
        print(col, nums, p_value)
        plt.figure(figsize=(12,10)) 
        
        for j in max_tod.index:
            plt.errorbar(decayBmal1_filtered[j], max_tod[j], xerr=decayBmal1_e[j], fmt='o', markersize=10, label=celllines[j])
            
        popt, pcov = curve_fit(linear_func, decayBmal1_filtered, max_tod)
        plt.plot(decayBmal1_filtered, linear_func(decayBmal1_filtered, *popt), '-k', label=str(col)+'\nCorr.={:.2f}, p-value={:.3f}'.format(corr_coef, p_value))
    
        plt.xlabel('Decay Bmal1 [1/hours]')
        plt.ylabel('Maximum relative ToD response')
        plt.legend(loc='best')
        # plt.savefig(path2+'Decay_max_ToD_'+str(col)+'.svg')
        plt.show()


