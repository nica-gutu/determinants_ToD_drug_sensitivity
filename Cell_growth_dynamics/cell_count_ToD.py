#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 10:08:58 2024

@author: nicagutu
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = '/Users/nicagutu/Nextcloud/Manuscripts/ToD_Confounders/Manuscript/to_upload/Data_toupload/Fig4_5_6/'
data = pd.read_excel(path+'271122_ToD_Complete.xlsx', sheet_name='Cisplatin')

time = data['Time']
averages = {}
errors = {}

ratios = []
column_names = []

for column in data.columns:
    if column.startswith('D'): ### 5FU: B, Cisplatin: D, Doxorubicin: F
        suffix = column[1:]  
        c_column = 'E' + suffix  ###  5FU: C, Cisplatin: E, Doxorubicin: G
        
        if c_column in data.columns:
            mean_values = (data[column] + data[c_column]) / 2
            std_errors = np.std([data[column], data[c_column]], axis=0)
            
            averages[suffix] = mean_values
            errors[suffix] = std_errors
            
            initial_value = mean_values.iloc[0]
            final_value = mean_values.iloc[-1]
            ratio = final_value/initial_value
            ratios.append(ratio)
            column_names.append(suffix)


plt.figure(figsize=(10, 6))
for suffix in averages:
    plt.errorbar(time, averages[suffix], yerr=errors[suffix], label=f'Dose {suffix}', capsize=5)
plt.xlabel('Time')
plt.ylabel('Cell counts')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(column_names, ratios, '-o', color='skyblue')
plt.xlabel('Column Suffix')
plt.ylabel('Final/Initial Ratio')
plt.show()


