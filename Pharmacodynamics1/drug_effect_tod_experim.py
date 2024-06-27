#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 16:39:37 2024

@author: nicagutu
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

def linear_func(x, m, b):
    return m * x + b

path ='Raw_Data/'
path2 = 'output/'
file = 'Circadian_properties_lumicycle'

data = pd.read_excel(path + file + '.xlsx', header=0, sheet_name='Maximum_ToD_Einf_h')
print(data)

#Heatmap with correlation coefficients
correlation_matrix = data.corr()
plt.figure(figsize=(22, 20))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", annot_kws={"size": 15},
            cbar_kws={"label": "Correlation Coefficient"})
plt.show()

celllines = data['Cellline']


#%%GRinf

col_list = ['Alpelisib', 'Paclitaxel', 'Alisertib']

plt.figure(figsize=(12, 10))

GRinf_all = []
max_tod_all = []

count = 0

markers = ['o', '*', 'D']

for col in col_list:
    GRinf = data[str(col) + '_Einf'].dropna()
    max_tod = data[col][GRinf.index]
    max_tod = (max_tod-min(max_tod))/(max(max_tod)-min(max_tod))
    GRinf = (GRinf-min(GRinf))/(max(GRinf)-min(GRinf))
    
    GRinf_all.append(GRinf.values)
    max_tod_all.append(max_tod.values)
    
    for j in max_tod.index:
        plt.plot(GRinf[j], max_tod[j], markers[count], markersize=10, label=f'{col} {celllines[j]}')
        
    count += 1
    
GRinf_all = np.concatenate(GRinf_all)
max_tod_all = np.concatenate(max_tod_all)

popt, pcov = curve_fit(linear_func, GRinf_all, max_tod_all)
corr_coef, p_value = pearsonr(GRinf_all, max_tod_all)

plt.plot(GRinf_all, linear_func(GRinf_all, *popt), '-k', label=f'Corr.={corr_coef:.2f}, p-value={p_value:.4f}')

plt.xlabel('Normalized GRinf')
plt.ylabel('Normalized maximum relative ToD response')
plt.legend(loc='best')
# plt.savefig(path2+'GRinf_max_ToD_more_drugs.svg')
plt.show()

#%%Hill coefficient

col_list = ['Adavosertib', 'Paclitaxel']

plt.figure(figsize=(12, 10))

hh_all = []
max_tod_all = []

markers = ['1', '*']
count = 0

for col in col_list:
    hh = data[str(col)+'_h'].dropna()
    max_tod = data[col][hh.index]
        
    max_tod = max_tod.dropna()
    hh = hh[max_tod.index]
    
    max_tod = (max_tod-min(max_tod))/(max(max_tod)-min(max_tod))
    hh = (hh-min(hh))/(max(hh)-min(hh))
    
    hh_all.append(hh.values)
    max_tod_all.append(max_tod.values)

    for j in max_tod.index:
        plt.plot(hh[j], max_tod[j], markers[count], markersize=10, label=f'{col} {celllines[j]}')
   
    count += 1
    
hh_all = np.concatenate(hh_all)
max_tod_all = np.concatenate(max_tod_all)

popt, pcov = curve_fit(linear_func, hh_all, max_tod_all)
corr_coef, p_value = pearsonr(hh_all, max_tod_all)

plt.plot(hh_all, linear_func(hh_all, *popt), '-k', label=f'Corr.={corr_coef:.2f}, p-value={p_value:.4f}')

plt.xlabel('Normalized Hill coefficient')
plt.ylabel('Normalized maximum relative ToD response')
plt.legend(loc='best')
# plt.savefig(path2+'Hillcoeff_max_ToD_more_drugs.svg')
plt.show()



