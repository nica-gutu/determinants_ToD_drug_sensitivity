#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:41:01 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod(A, t, T, phi): #circadian modulation
    return Amp*np.sin(t*(2*np.pi/T)+phi)

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

path ='output/'

#Simulation parameters
tf = 120
t = np.linspace(0, tf, tf)
toD = np.arange(0, 25, step=1) #time of day
step = 24
toE = np.arange(step, tf+step, step=step)

colors = plt.cm.tab20(range(len(toE)))

#Growth parameters
division_time = 24
k = np.log(2)/division_time #growth rate

#Drug parameters
Sm = 1 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect
x_0 = 1

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

#%%No circadian modulation
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c = 10 ** exp_factors

colors2 = plt.cm.tab20b(np.linspace(0, 1, len(c)))
fig = plt.figure(figsize = (10,10))
for j in range(len(c)):
    plt.plot(t, growthcurve(x_0, t, k, kl, LC50, Sm, c[j], SC50, h), color=colors2[j])
    if j%5==0:
        plt.plot(t, growthcurve(x_0, t, k, kl, LC50, Sm, c[j], SC50, h), color=colors2[j], label='{:.3f}'.format(c[j]))
for i in range(len(toE)):
    plt.axvline(x=toE[i], color='grey', linestyle='--')  
plt.xlabel('Time')
plt.ylabel('Cell count')
plt.legend(loc='upper left', title='Concentration')
# plt.savefig(path+'Cell_count_di1v_time'+str(division_time)+'.svg')
plt.show()

hGR_values = []
inf_values = []
EC50_values = []
# hGR_slopes = []
fig = plt.figure(figsize = (10,10))
for i in range(len(toE)):
    drc = []
    for j in range(len(c)):
        drc.append(growthcurve(x_0, toE[i], k, kl, LC50, Sm, c[j], SC50, h)/growthcurve(x_0, toE[i], k, kl, LC50, Sm, 0, SC50, h))    
    plt.plot(c, drc, color=colors[i], label=str(toE[i])+'h')
    popt, pcov = curve_fit(GR, c, drc)#, bounds=([0.9,0,1], [1.1,0.5,3]))
    inf_values.append(popt[0])
    EC50_values.append(popt[1])
    hGR_values.append(popt[2])
    plt.plot(popt[1], growthcurve(x_0, toE[i], k, kl, LC50, Sm, popt[1], SC50, h)/growthcurve(x_0, toE[i], k, kl, LC50, Sm, 0, SC50, h), 'o', color=colors[i])
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Normalized cell count to untreated')
plt.legend(loc='best', title='Time of eval')
# plt.savefig(path+'Norm_cell_count_untr_div_time'+str(division_time)+'.svg')
plt.show()


fig = plt.figure(figsize = (10,10))
for i in range(len(toE)):
    plt.plot(toE[i], hGR_values[i], 'o', color=colors[i], label=str(toE[i])+'h')
plt.xlabel('Time of Evaluation')
plt.ylabel('Hill coefficient')
plt.legend(loc='best', title='Time of eval')
# plt.savefig(path+'hGR_vs_ToE_div_time'+str(division_time)+'.svg')
plt.show()


fig = plt.figure(figsize = (10,10))
for i in range(len(toE)):
    plt.plot(toE[i], inf_values[i], 'o', color=colors[i], label=str(toE[i])+'h')
plt.xlabel('Time of Evaluation')
plt.ylabel('GRinf')
plt.legend(loc='best', title='Time of eval')
# plt.savefig(path+'GRinf_vs_ToE_div_time'+str(division_time)+'.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
for i in range(len(toE)):
    plt.plot(toE[i], EC50_values[i], 'o', color=colors[i], label=str(toE[i])+'h')
plt.xlabel('Time of Evaluation')
plt.ylabel('GR50')
plt.legend(loc='best', title='Time of eval')
# plt.savefig(path+'GREC50_vs_ToE_div_time'+str(division_time)+'.svg')
plt.show()

#%%Concentration with circadian modulation
toD = np.linspace(0, 24, 200) #time of day

division_time = 24
k = np.log(2)/division_time #growth rate

#Drug parameters
Sm = 0 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect
x_0 = 1

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

mr = []
fig = plt.figure(figsize = (10,10))
for i in range(len(toE)):
    drc = []
    EC50 = EC50_values[-1]
    h = hGR_values[i]
    Sm = inf_values[-1]
    c = EC50*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c)):
        drc.append(growthcurve(x_0, toE[i], k, kl, LC50, Sm, c[j], SC50, h))
    drc = drc/drc[0]
    plt.plot(toD, drc, color=colors[i], linewidth=5, label=str(toE[i])+'h')
    mr.append(abs(max(drc)-min(drc)))
plt.xticks([0,6,12,18,24])
plt.xlabel('Time [hours]')
plt.ylabel('ToD response')
plt.legend(loc='best', title='Time of eval')
# plt.savefig(path+'ToD_ToE_div_time'+str(division_time)+'.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
plt.plot(toE, mr)
plt.xlabel('Time [hours]')
plt.ylabel('Maximum range of ToD')
# plt.savefig(path+'Max_range_ToD_mod_concentr_diff_ToE_div_time'+str(division_time)+'.svg')
plt.show()
