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

path ='output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

#Simulation parameters
tf = 145
t = np.linspace(0, tf, tf)
toD = np.arange(0, 24.1, step=0.1) #time of day
colors2 = plt.cm.tab20b(range(len(toD)))   
 
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate

#Drug parameters
Sm = 0 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0.05 #killing effect
LC50 = 1 #half killing effect

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

opt = 'change_c'
exp_factors = np.array(np.arange(-4, 2.2, step=0.2), dtype=float) 
c1 = (10**exp_factors) 

#Find initial EC50
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

#%%Survival curve
exp_factors = np.array(np.arange(-4, 2.2, step=0.2), dtype=float) 
c1 = (10**exp_factors) 

# var_vect = [-ini, -ini+0.0005, -ini+0.0025, -ini+0.01, -ini+0.03, 0, 0.02, 0.05, 0.2, 1, 3, 10]
var_vect = np.arange(-ini, 10, step=0.01)
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
    c2 = (popt[1]+var_vect[ii])*circmod_c(modul_strgth, Amp, toD, T, phi)
    conc.append(popt[1]+var_vect[ii])
    for j in range(12,13):
        drc = (growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
        plt.plot(c2[j], drc, 'o', markersize=15, color=colors[ii], label='{:.2f}'.format(popt[1]+var_vect[ii]))
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Survival fraction')
plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Survival_modulated_concentration_'+str(opt)+'.svg')
plt.show()

#%%Effective concentration
fig = plt.figure(figsize = (10,10))
for j in range(len(var_vect)):
    EC50 = EC50_values[j]
    c = (EC50+var_vect[j])*circmod_c(modul_strgth, Amp, toD, T, phi)
    plt.plot(toD, c, markersize=15, label='{:.3f}'.format(var_vect[j]))
plt.xticks([0,6,12,18,24])
plt.xlabel('Time [hrs]')
plt.ylabel('Modulated Concentration')
plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Modul_concentration_'+str(opt)+'.svg')
plt.show()

#%%ToD response
mr = []
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    ratio = []
    c2 = (EC50_values[ii]+var_vect[ii])*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c2)):
        ratio.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    plt.plot(toD, ratio/ratio[0], linewidth=4, color=colors[ii], label='{:.3f}'.format(EC50_values[ii]+var_vect[ii])) ###change label
    mr.append(abs(max(ratio/ratio[0])-min(ratio/ratio[0])))
plt.xticks([0,6,12,18,24])
plt.ylabel('ToD response')
plt.xlabel('Time of day [hours]')
plt.legend(loc='best', title='Stress strength')
# plt.savefig(path+'Timeofday_'+str(opt)+'.svg')
plt.show()


#%%Maximum range of ToD response
plt.figure(figsize = (12,10))
plt.plot(conc, mr, 'o', markersize=15, linewidth=5)
plt.xscale('log')
plt.ylabel('Maximum range ToD response')
plt.xlabel('Concentration')  
# plt.savefig(path+'Max_ToD_'+str(opt)+'.svg')
plt.show()

