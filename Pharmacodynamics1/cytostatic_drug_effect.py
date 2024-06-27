#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:17:16 2023

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
Sm = 1 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect

opt = 'change_h' ####change label
var_vect =  [1, 1.5, 2, 3, 5]#, 10] #Changing Hill coefif
# var_vect = [0.5, 1, 1.5, 2, 2.5, 5, 10] #Changing Sm
###### var_vect = [0.5, 1, 1.5, 2, 5, 10] #Changing SC50 #####no significant changes
colors = plt.cm.tab20(range(len(var_vect)))

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

#%%Modulated concentration and cell count
# c = 1*circmod_c(modul_strgth, Amp, toD, T, phi)
# fig = plt.figure(figsize = (10,10))
# for j in range(len(toD)):
#     plt.plot(toD[j], c[j], 'o', markersize=15, color=colors2[j])
# plt.xlabel('Time [hrs]')
# plt.ylabel('Modulated Concentration')
# plt.title('Hill coeff:'+str(h)+' kl:'+str(kl)+' LC50:'+str(LC50))
# plt.savefig(path+'Modul_concentration_'+str(opt)+'.svg')
# plt.show()

# fig = plt.figure(figsize = (10,10))
# for j in range(len(c)):
#     # Sm = 5
#     # SC50 = 5
#     h = 5
#     plt.plot(t, growthcurve(1, t, k, kl, LC50, Sm, c[j], SC50, h), linewidth=5, color=colors2[j])
# plt.xlabel('Time [hrs]')
# plt.ylabel('Normalized cell count')
# plt.title('Hill coeff:'+str(h)+' Sm:'+str(Sm)+' SC50:'+str(SC50))
# # plt.savefig(path+'Cellcount_'+str(opt)+str(h)+'_doublingtime_'+str(doubling_time1)+'.svg') #####change
# plt.show()

#%%Survival curve
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

plt.figure(figsize = (12,10))
EC50_values = []
for ii in range(len(var_vect)):
    drc = []
    h = var_vect[ii] ######change variable
    for j in range(len(c1)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    plt.plot(c1, drc, linewidth=4, color=colors[ii], label='h:'+str(var_vect[ii])) ######change label
    popt, pcov = curve_fit(GR, c1, drc)
    EC50_values.append(popt[1])
    c2 = popt[1]*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c2)):
        drc = (growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
        plt.plot(c2[j], drc, 'o', markersize=15, color=colors2[j])
plt.xscale(value='log')
plt.xlabel('Concentration')
plt.ylabel('Dose response curve')
plt.legend(loc='best')
# plt.title('Changing Sm; Hill coeff:'+str(h)+' SC50:'+str(SC50)) ####comment this line or the next one
# plt.title('Changing h; Sm:'+str(Sm)+' SC50:'+str(SC50)) ####comment
# plt.savefig(path+'DRC_modulated_concentration_'+str(opt)+'_doublingtime_'+str(doubling_time1)+'.svg')
plt.show()

#%%ToD response
plt.figure(figsize = (12,10))
for ii in range(len(var_vect)):
    ratio = []
    h = var_vect[ii] ######change variable
    c2 = EC50_values[ii]*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c2)):
        ratio.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    plt.plot(toD, ratio/ratio[0], linewidth=4, color=colors[ii], label='h:'+str(var_vect[ii])) #####change label
plt.ylabel('Drug sensitivity')
plt.xlabel('Time of day [hours]')
plt.legend(loc='best')
# plt.title('Changing Sm; Hill coeff:'+str(h)+' SC50:'+str(SC50)) ####comment this line or the next one
# plt.title('Changing h; Sm:'+str(Sm)+' SC50:'+str(SC50)) ####comment
# plt.savefig(path+'Timeofday_'+str(opt)+'_doublingtime_'+str(doubling_time1)+'.svg')
plt.show()



