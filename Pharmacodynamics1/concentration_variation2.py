#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:27:16 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings("ignore")

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

path ='Output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return modul_strgth*(Amp*np.sin(t*(2*np.pi/T)+phi))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

#Simulation parameters
tf = 145
t = np.linspace(0, tf, tf)
toD = np.arange(0, 28, step=4) #time of day
colors2 = plt.cm.tab20b(range(len(toD)))   
 
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate
Sm = 0 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0.05 #killing effect
LC50 = 1 #half killing effect

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.1

opt = 'change_c'
exp_factors = np.array(np.arange(-4, 2.2, step=0.2), dtype=float) 
c1 = (10**exp_factors) 

plt.figure(figsize = (12,10))
ini = 0
drc = []
for j in range(len(c1)):
    drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
plt.plot(c1, drc, linewidth=4)#, label='{:.3f}'.format(var_vect[ii])) ###change label
popt, pcov = curve_fit(GR, c1, drc)
ini = popt[1]
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Survival fraction')
plt.show()
print(ini)


exp_factors = np.array(np.arange(-4, 2.2, step=0.1), dtype=float) 
c1 = (10**exp_factors) 

var_vect = [0.001, 0.01, 0.1, 1, 10]
var_vect = c1
colors = plt.cm.tab20(range(len(var_vect)))

plt.figure(figsize = (12,10))
EC50_values = []
E_inf  = []
conc = []
for ii in range(len(var_vect)):
    drc = []
    for j in range(len(c1)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c1, drc, linewidth=4, color=colors[ii])#, label='{:.3f}'.format(var_vect[ii])) ###change label
    popt, pcov = curve_fit(GR, c1, drc)
    EC50_values.append(popt[1])
    E_inf.append(popt[0])
    c2 = (var_vect[ii])*10**(circmod_c(modul_strgth, Amp, toD, T, phi))
    conc.append(var_vect[ii])
    for j in range(len(toD)):
        drc = (growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
        plt.plot(c2[j], drc, 'o', markersize=15, color=colors[ii])#, label='{:.3f}'.format(var_vect[ii]))
    plt.plot(c2[0], drc, 'o', markersize=15, color=colors[ii], label='{:.4f}'.format(var_vect[ii]))
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Survival fraction')
# plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Survival_modulated_concentration_'+str(opt)+'.svg')
plt.show()

toD = np.arange(0, 24.1, step=.1) #time of day
fig = plt.figure(figsize = (10,10))
for j in range(len(var_vect)):
    EC50 = EC50_values[j]
    c = (var_vect[j])*10**(circmod_c(modul_strgth, Amp, toD, T, phi))
    plt.plot(toD, c, markersize=15, label='{:.3f}'.format(var_vect[j]))
plt.xticks([0,6,12,18,24])
plt.xlabel('Time [hrs]')
plt.ylabel('Modulated Concentration')
# plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Modul_concentration_'+str(opt)+'.svg')
plt.show()


mr = []
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    ratio = []
    c2 = (var_vect[ii])*10**circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c2)):
        ratio.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    plt.plot(toD, (ratio-np.mean(ratio)+5)/5, linewidth=4, color=colors[ii], label='{:.3f}'.format(var_vect[ii])) ###change label
    mr.append(abs(max(ratio)-min(ratio)))
plt.xticks([0,6,12,18,24])
plt.ylabel('Relative ToD response')
plt.xlabel('Time of day [hours]')
# plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Timeofday_'+str(opt)+'.svg')
plt.show()

plt.figure(figsize = (12,10))
plt.plot(var_vect, mr, 'o', markersize=15, linewidth=5)
plt.xscale('log')
plt.ylabel('Maximum range ToD response')
plt.xlabel('Concentration')  #######change label
#plt.savefig(path+'Max_ToD_'+str(opt)+'.svg')
plt.show()

