#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 13:10:29 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path ='output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h)))#-t*(kl*c**h)/(LC50**h+c**h))

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)


#Simulation parameters
tf = 145
t = np.linspace(0, tf, tf)
toD = np.arange(0, 26, step=2) #time of day
    
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate
step = 4
doubling_times = np.arange(12, 48+step, step=step)

#Drug parameters
Sm = 1 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

c = 1*circmod_c(modul_strgth, Amp, toD, T, phi)
x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
# print(x0_vect)
colors = plt.cm.tab20b(range(len(c)))

fig = plt.figure(figsize = (10,10))
for j in range(len(toD)):
    plt.plot(toD[j], c[j], 'o', markersize=15, color=colors[j])
plt.xlabel('Time [hours]')
plt.ylabel('Modulated Concentration')
# plt.savefig(path+'Modul_concentration.svg')
plt.show()

#%%Normalized cell count
fig = plt.figure(figsize = (10,10))
for j in range(len(x0_vect)):
    # time_adj = np.arange(0, tf-toD[j], step=1)
    plt.plot(t, growthcurve(x0_vect[j], t, k, kl, LC50, Sm, c[j], SC50, h), linewidth=5, color=colors[j])
plt.xlabel('Time [hours]')
plt.ylabel('Normalized cell count')
# plt.savefig(path+'Cellcount.svg')
plt.show()

#%%Survival fraction
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

fig = plt.figure(figsize = (10,10))
drc = []
for j in range(len(c1)):
    drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
plt.plot(c1, drc, linewidth=4)
popt, pcov = curve_fit(GR, c1, drc)
EC50_value = popt[1]
c2 = EC50_value*circmod_c(modul_strgth, Amp, toD, T, phi)
for j in range(len(c2)):
    drc = (growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c2[j], drc, 'o', markersize=15, color=colors[j])
plt.xscale(value='log')
# plt.xlim([0,3])
plt.xlabel('Concentration')
plt.ylabel('Survival fraction')
plt.legend(loc='best')
# plt.savefig(path+'DRC_modulated_concentration.svg')
plt.show()

#%%ToD response
fig = plt.figure(figsize = (10,10))
for j in range(len(c2)):
    norm = growthcurve(1, t, k, kl, LC50, Sm, c2[0], SC50, h)[-1]
    ratio = growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]
    plt.plot(toD[j], ratio/norm, 'o', markersize=15, color=colors[j])
plt.ylabel('ToD response')
plt.xlabel('Time of day [hrs]')
# plt.savefig(path+'Timeofday.svg')
plt.show()


