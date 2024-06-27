#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:13:07 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

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
t = np.linspace(0, tf, tf) #entire time of the simulation
toD = np.arange(0, 24.01, step=0.01) #time of day
    
#Growth parameters
doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate

#Drug parameters
Sm = 0. #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0.01 #killing effect
LC50 = 1 #half killing effect

#Simulation range
opt = 'change_T' #####change label
# var_vect = [0.1, 0.2, 0.5, 0.75, 1] #Changing Amp #np.arange(0, 1.025, step=0.025)#
var_vect =  [16, 24, 32, 40, 48] #Changing T #np.arange(16, 48, step=1)#

colors = plt.cm.tab20(range(len(var_vect)))

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 1

#%%Effective concentration
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    T = var_vect[ii] ####change variable
    c = 1*circmod_c(modul_strgth, Amp, toD, T, phi)        
    plt.plot(toD, c, linewidth = 4, color=colors[ii], label='A:'+str(var_vect[ii])) ######change label
plt.xlabel('Time of day [hours]')
plt.xticks([0,6,12,18,24])
plt.ylabel('Modulated Concentration')
plt.legend(loc='best')
# plt.savefig(path+'Modul_concentration_'+str(opt)+'.svg')
plt.show()

#%%Survival curve
exp_factors = np.array(np.arange(-7, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

plt.figure(figsize = (12,10))
EC50_values = []
for ii in range(len(var_vect)):
    drc = []
    T = var_vect[ii] ######change variable
    for j in range(len(c1)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c1, drc, linewidth=4, color=colors[ii], label='T:'+str(var_vect[ii])) ######change label
    popt, pcov = curve_fit(GR, c1, drc)
    EC50_values.append(popt[1])
    c2 = popt[1]*circmod_c(modul_strgth, Amp, toD, T, phi)
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
mr = []
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    ratio = []
    T = var_vect[ii] ######change variable
    c2 = EC50_values[ii]*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c2)):
        ratio.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    plt.plot(toD, ratio/ratio[0], linewidth=4, color=colors[ii], label='T:'+str(var_vect[ii])) #####change label
    mr.append(abs(max(ratio/ratio[0])-min(ratio/ratio[0])))
plt.xticks([0,6,12,18,24])
plt.ylabel('ToD response')
plt.xlabel('Time of day [hours]')
plt.legend(loc='best')
# plt.savefig(path+'Timeofday_'+str(opt)+'_doublingtime_'+str(doubling_time1)+'.svg')
plt.show()

#%%Maximum range of ToD response
plt.figure(figsize = (12,10))
plt.plot(var_vect, mr, '-o', markersize=15, linewidth=5)
corr_coef, p_value = pearsonr(var_vect, mr)
print(corr_coef)
plt.ylabel('Maximum range relative ToD response')
plt.xlabel('Period [hours]') #######Change
# plt.savefig(path+'Max_ToD_'+str(opt)+'_doublingtime_'+str(doubling_time1)+'.svg')
plt.show()



