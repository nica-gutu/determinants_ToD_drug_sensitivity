#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:36:27 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

path ='output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def pulse(time, b):
    t_normalized = time % (2 * b)
    result = np.where(t_normalized < b, 0, 1)
    return result
    
def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

#Simulation parameters
tf = 60
t = np.linspace(0, tf, tf)
toD = np.arange(0, 25, step=1) #time of day

colors2 = plt.cm.tab20b(range(len(toD)))
    
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate

#Drug parameters
Sm = 0. #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0.01 #killing effect
LC50 = 1 #half killing effect

opt = 'change_b'
var_vect =  [6,12]
colors = plt.cm.tab20(range(len(var_vect)))

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 1

#%%Pulse function 
plt.figure(figsize=(10,8))
plt.plot(t, pulse(t, 12), '-', linewidth=5)
plt.xlabel('Time')
plt.ylabel('f(t)')
# plt.savefig(path+'Pulse_function_sketch.svg')
plt.show()

#%%Modulated concentration
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    b = var_vect[ii]
    c = 1*pulse(toD, b)     
    plt.plot(toD, c, linewidth=5, label='b:'+str(var_vect[ii])) ######change
plt.xlabel('Time of day [hours]')
plt.ylabel('Modulated Concentration')
plt.legend(loc='best')
# plt.savefig(path+'Modul_concentration_pulse_function.svg')
plt.show()

#%%Survival curve
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

plt.figure(figsize = (12,10))
EC50_values = []
for ii in range(len(var_vect)):
    drc = []
    b = var_vect[ii] 
    for j in range(len(c1)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c1, drc, linewidth=4, color=colors[ii], label='b:'+str(var_vect[ii])) 
    popt, pcov = curve_fit(GR, c1, drc)
    EC50_values.append(popt[1])
    c2 = popt[1]*pulse(toD, b)
    drc = []
    for j in range(len(c2)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c2, drc, 'o', markersize=15, color=colors[ii])
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Survival fraction')
plt.legend(loc='best')
# plt.savefig(path+'DRC_modulated_concentration_'+str(opt)+'_doublingtime_'+str(doubling_time1)+'.svg')
plt.show()

#%%ToD response
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    ratio = []
    b = var_vect[ii] 
    c2 = EC50_values[ii]*pulse(toD, b)
    for j in range(len(c2)):
        ratio.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    plt.plot(toD, ratio/ratio[0], linewidth=5, label='b:'+str(var_vect[ii])) 
plt.ylabel('Drug sensitivity')
plt.xlabel('Time of day [hours]')
plt.legend(loc='best')
# plt.savefig(path+'Timeofday_pulse_function.svg')
plt.show()
