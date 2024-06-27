#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:04:50 2023

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
toD = np.arange(0, 25, step=1) #time of day
    
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

colors = plt.cm.tab20(range(len(doubling_times)))

#%%No circadian modulation
c = 1
x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)

#Adjusting time
fig = plt.figure(figsize = (10,10))
plt.plot(t, growthcurve(1, t, k, kl, LC50, Sm, c, SC50, h), linewidth=5)
for j in range(1,len(x0_vect)):
    time_adj = np.arange(toD[j], tf, step=1)
    plt.plot(time_adj[0], x0_vect[j], 'o', markersize=9, label=str(toD[j])+'h')
plt.xlabel('Time')
plt.ylabel('Normalized cell count')
plt.title('Division time = '+str(doubling_time1)+' h')
plt.legend(loc='best', title='Treatment time')
# plt.savefig(path+'Cell_count_division'+str(doubling_time1)+'fixed_conc.svg')
plt.show()


fig = plt.figure(figsize = (10,10))
for j in range(len(x0_vect)):
    time_adj = np.arange(0, tf-toD[j], step=1)
    plt.plot(time_adj, (growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c, SC50, h)), linewidth=5)
plt.xlabel('Time [hrs]')
plt.ylabel('Normalized cell count')
plt.title('Division time = '+str(doubling_time1)+' h')
# plt.savefig(path+'Cell_count_division'+str(doubling_time1)+'early_start_fixed_conc.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
for j in range(len(x0_vect)):
    time_adj = np.arange(0, tf-toD[j], step=1)
    plt.plot(time_adj, (growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c, SC50, h)), linewidth=5)
plt.xlabel('Time [hrs]')
plt.ylabel('Normalized cell count')
plt.title('Division time = '+str(doubling_time1)+' h')
plt.xlim([-5,120])
# plt.savefig(path+'Cell_count_division'+str(doubling_time1)+'early_start_fixed_conc.svg')
plt.show()

#Wrong decision
fig = plt.figure(figsize = (10,10))
for i in range(len(doubling_times)):
    k = np.log(2)/doubling_times[i]
    drc = []
    x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
    for j in range(len(x0_vect)):
        time_adj = np.arange(toD[j], tf+1, step=1)
        drc.append(growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c, SC50, h)[-1])#/growthcurve(x0_vect[j], t, k, kl, LC50, Sm, c, SC50, h)[0])
    plt.plot(toD, drc/drc[0], color=colors[i], label=str(doubling_times[i])+'h', linewidth=5)
# plt.ylim([1,2])
# plt.xlim([0,24])
plt.xlabel('Time of treatment')
plt.ylabel('Final relative response')
plt.legend(loc='best', title='Division times')
# plt.savefig(path+'ToD_wrong_decision_fixed_concentration.svg')
plt.show()

#Solution
fig = plt.figure(figsize = (10,10))
for i in range(len(doubling_times)):
    k = np.log(2)/doubling_times[i]
    drc = []
    x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
    for j in range(len(x0_vect)):
        time_adj = np.arange(toD[j], tf+1, step=1)
        x0 = growthcurve(x0_vect[j], t, k, kl, LC50, Sm, c, SC50, h)[toD[j]]
        drc.append(growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c, SC50, h)[-1]/x0_vect[j])
    plt.plot(toD, drc/drc[0], color=colors[i], label=str(doubling_times[i])+'h', linewidth=5) 
# plt.ylim([0,1])
plt.xlim([0,24])
plt.xlabel('Time of treatment')
plt.ylabel('Final relative response')
plt.legend(loc='best', title='Division times')
# plt.savefig(path+'ToD_right_decision_fixed_concentration.svg')
plt.show()


#%%Circadian modulated concentration

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5


toD = np.arange(0, 25, step=1) #time of day
c = 1*circmod_c(modul_strgth, Amp, toD, T, phi)
x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
# print(x0_vect)
colors2 = plt.cm.tab20b(range(len(toD)))

#Modulated concentration
fig = plt.figure(figsize = (10,10))
plt.scatter(toD, c, c=colors2, s=100)
plt.xlabel('Time [hours]')
plt.ylabel('Modulated Concentration')
# plt.savefig(path+'Modulated_concentration.svg')
plt.show()

doubling_time1 = 24
k = np.log(2)/doubling_time1 #growth rate

#Adjusting time
fig = plt.figure(figsize = (10,10))
plt.plot(t, growthcurve(1, t, k, kl, LC50, Sm, 1, SC50, h))
for j in range(len(x0_vect)):
    time_adj = np.arange(toD[j], tf, step=1)
    plt.plot(time_adj, growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c[j], SC50, h), label=str(toD[j])+'h')
plt.xlabel('Time [hours]')
plt.ylabel('Cell count')
plt.title('Division time = '+str(doubling_time1)+' h')
# plt.legend(loc='best', title='Treatment time')
plt.savefig(path+'Cell_count_division'+str(doubling_time1)+'late_start_mod_concentr.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
for j in range(len(x0_vect)):
    time_adj = np.arange(0, tf-toD[j], step=1)
    plt.plot(time_adj, growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c[j], SC50, h), linewidth=5)
plt.xlabel('Time [hours]')
plt.ylabel('Normalized cell count')
plt.title('Division time = '+str(doubling_time1)+' h')
# plt.savefig(path+'Cell_count_division'+str(doubling_time1)+'early_start_mod_concentr.svg')
plt.show()

#Survival curve
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

EC50_values = []
fig = plt.figure(figsize = (10,10))
for i in range(len(doubling_times)):
    drc = []
    k = np.log(2)/doubling_times[i]
    for j in range(len(c1)):
        drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
    popt, pcov = curve_fit(GR, c1, drc)
    EC50_values.append(popt[1])
    plt.plot(c1, drc, linewidth=5, label=str(doubling_times[i])+'h')
    # plt.plot(popt[1], growthcurve(1, t, k, kl, LC50, Sm, popt[1], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1], 'o')
    c = popt[1]*circmod_c(modul_strgth, Amp, toD, T, phi)
    for j in range(len(c)):
        drc2 = (growthcurve(1, t, k, kl, LC50, Sm, c[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
        plt.plot(c[j], drc2, 'o', markersize=10)#, color=colors2[j])
plt.xscale(value='log')
plt.ylabel('Survival')
plt.xlabel('Dose ')
plt.legend(loc='best', title='Division times')
# plt.savefig(path+'DRC_EC50_modulated_concentration.svg')
plt.show()

#Wrong decision
fig = plt.figure(figsize = (10,10))
for i in range(len(doubling_times)):
    k = np.log(2)/doubling_times[i]
    drc = []
    c = EC50_values[i]*circmod_c(modul_strgth, Amp, toD, T, phi)
    x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
    for j in range(len(x0_vect)):
        time_adj = np.arange(toD[j], tf+1, step=1)
        drc.append(growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c[j], SC50, h)[-1])
    plt.plot(toD, drc/drc[0], color=colors[i], linewidth=5, label=str(doubling_times[i])+'h') 
plt.xlabel('Time of day (~initial density)')
plt.ylabel('Final relative response to t=0')
plt.legend(loc='best', title='Division times')
# plt.savefig(path+'ToD_wrong_decision_modulated_concentration.svg')
plt.show()

#Solution 
fig = plt.figure(figsize = (10,10))
for i in range(len(doubling_times)):
    k = np.log(2)/doubling_times[i]
    c = EC50_values[i]*circmod_c(modul_strgth, Amp, toD, T, phi)
    x0_vect = growthcurve(1, toD, k, kl, LC50, Sm, c, SC50, h)
    drc = []
    for j in range(len(x0_vect)):
        time_adj = np.arange(toD[j], tf+1, step=1)
        drc.append(growthcurve(x0_vect[j], time_adj, k, kl, LC50, Sm, c[j], SC50, h)[-1]/x0_vect[j])
    plt.plot(toD, drc/drc[0], color=colors[i], linewidth=5, label=str(doubling_times[i])+'h')
plt.xlabel('Time of day (~initial density)')
plt.ylabel('Final relative response to initial density')
plt.legend(loc='best', title='Division times')
# plt.savefig(path+'ToD_right_decision_modulated_concentration.svg')
plt.show()




