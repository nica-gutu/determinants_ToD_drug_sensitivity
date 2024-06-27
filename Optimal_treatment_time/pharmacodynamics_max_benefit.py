#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 10:03:14 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'

path ='Output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return (1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

def c_exp(kdecay, t): #exponential decay (drug stability)
    return np.exp(-kdecay*t)

#Simulation parameters
tf = 145
t = np.linspace(0, tf, tf)
toD = np.arange(0, 24.01, step=0.01) #time of day

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

##Compute IC50
exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
c1 = 10 ** exp_factors

drc = []
for j in range(len(c1)):
    drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
popt, pcov = curve_fit(GR, c1, drc)
IC50 = popt[1]

#Circadian parameters healthy tissue
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5
c2 = IC50*circmod_c(modul_strgth, Amp, toD, T, phi) ##no damping

tod_healthy = []
for j in range(len(c2)):
    tod_healthy.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
tod_healthy = tod_healthy/tod_healthy[0]

time_max_health = np.argmax(tod_healthy)
time_min_health = np.argmin(tod_healthy)


plt.figure(figsize = (12,10))
plt.plot(toD, tod_healthy, linewidth=4, color='green', label='Healthy')

#%%Change Maximal drug effect
# #Circadian parameters tumor
# Amp = 1
# T = 24
# phi = 0 #0 or pi
# modul_strgth = 0.5
# d = 0

# #Parameters to change
# opt = 'change_kl'
# var_vect = np.arange(0, 0.01, step=0.00025)#[0.001, 0.0025, 0.005, 0.0075] #Changing kl # 
# colors = plt.cm.tab20c(range(len(var_vect)))

# exp_factors = np.array(np.arange(-3, 1.2, step=0.2), dtype=float) 
# c1 = 10 ** exp_factors

# EC50_values = []
# E_inf  = []
# for ii in range(len(var_vect)):
#     drc = []
#     kl = var_vect[ii] 
#     for j in range(len(c1)):
#         drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
#     popt, pcov = curve_fit(GR, c1, drc)
#     EC50_values.append(popt[1])
#     E_inf.append(popt[0])

# max_benefit = []
# min_benefit = []

# fold_change_max = []
# fold_change_min = []

# time_max = []
# time_min = []

# for ii in range(len(var_vect)):
#     tod_tumor = []
    
#     kl = var_vect[ii]

#     c2 = EC50_values[ii]*circmod_c(modul_strgth, Amp, toD, T, phi) 
#     for j in range(len(c2)):
#         tod_tumor.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
#     tod_tumor = tod_tumor/tod_tumor[0]
    
#     plt.plot(toD, tod_tumor, linewidth=4, color=colors[ii], label='GRinf: {:.2f}'.format(E_inf[ii])) #####change label
    
#     max_benefit.append((max(tod_healthy)-tod_tumor[time_max_health]))
#     min_benefit.append((min(tod_healthy)-tod_tumor[time_min_health]))

#     # fold_change_max.append((max(tod_healthy)/(tod_tumor[time_max_health])))
#     # fold_change_min.append((min(tod_healthy)/(tod_tumor[time_min_health])))

#     fold_change = tod_healthy/tod_tumor
#     time_max_benefit = toD[np.argmax(fold_change)]
#     time_min_benefit = toD[np.argmin(fold_change)]
#     fold_change_max.append(max(fold_change))
#     fold_change_min.append(min(fold_change))
#     time_max.append(time_max_benefit)
#     time_min.append(time_min_benefit)


# plt.xticks([0,6,12,18,24])
# plt.ylabel('ToD response')
# plt.xlabel('Time of day [hours]')
# # plt.legend(loc='best')
# # plt.savefig(path+'ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()


# # plt.figure(figsize = (12,10))
# # plt.plot(E_inf, max_benefit, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
# # plt.plot(E_inf, min_benefit, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
# # plt.ylabel('Benefit range of ToD response\n healthy tissue vs tumor')
# # plt.xlabel('GRinf') 
# # plt.legend(loc='best')
# # # plt.savefig(path+'Rel_max_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# # plt.show()

# plt.figure(figsize = (12,10))
# plt.plot(E_inf, fold_change_max, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
# plt.plot(E_inf, fold_change_min, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
# plt.ylabel('Fold-change of ToD response\n healthy tissue vs tumor')
# plt.xlabel('Emax') 
# plt.legend(loc='best')
# plt.savefig(path+'Fold-change_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111, polar=True)
# cmap = plt.get_cmap('YlGn')
# r = np.ones_like(time_max)
# colors = cmap(var_vect)
# scatter = ax.scatter(time_max, r, c=var_vect, s=200, cmap=cmap)
# ax.grid(False)  # Hide the grid
# ax.set_yticklabels([])  # Hide the radial labels
# hour_labels = [f'{i}h' for i in range(0, 24, 3)] 
# ax.set_xticks(np.linspace(0, 2 * np.pi, len(hour_labels), endpoint=False))
# ax.set_xticklabels(hour_labels)
# cbar = plt.colorbar(scatter, ax=ax, orientation='vertical')
# cbar.set_label('Emax')
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1)
# plt.title('Time of maximal benefit')
# plt.savefig(path+'Time_max_benefit_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()


#%%Change Hill coefficient
# #Circadian parameters tumor
# Amp = 1
# T = 24
# phi = 0 #0 or pi
# modul_strgth = 0.5
# d = 0

# #Parameters to change
# opt = 'change_h'
# var_vect = np.arange(0.1, 1, step=0.02)# [0.1, 0.25, 0.5, 0.75]#[1.5, 2, 3, 4, 5] #Changing h 
# colors = plt.cm.tab20c(range(len(var_vect)))

# exp_factors = np.array(np.arange(-3, 3.2, step=0.2), dtype=float) 
# c1 = 10 ** exp_factors

# EC50_values = []
# for ii in range(len(var_vect)):
#     drc = []
#     h = var_vect[ii] 
#     for j in range(len(c1)):
#         drc.append(growthcurve(1, t, k, kl, LC50, Sm, c1[j], SC50, h)[-1]/growthcurve(1, t, k, kl, LC50, Sm, 0, SC50, h)[-1])    
#     popt, pcov = curve_fit(GR, c1, drc)
#     EC50_values.append(popt[1])

# max_benefit = []
# min_benefit = []

# fold_change_max = []
# fold_change_min = []

# time_max = []
# time_min = []

# for ii in range(len(var_vect)):
#     tod_tumor = []
    
#     h = var_vect[ii]

#     c2 = EC50_values[ii]*circmod_c(modul_strgth, Amp, toD, T, phi) 
#     for j in range(len(c2)):
#         tod_tumor.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
#     tod_tumor = tod_tumor/tod_tumor[0]
    
#     plt.plot(toD, tod_tumor, linewidth=4, color=colors[ii], label='h:{:.2f}'.format(var_vect[ii])) #####change label
    
#     max_benefit.append((max(tod_healthy)-tod_tumor[time_max_health]))
#     min_benefit.append((min(tod_healthy)-tod_tumor[time_min_health]))

#     # fold_change_max.append((max(tod_healthy)/(tod_tumor[time_max_health])))
#     # fold_change_min.append((min(tod_healthy)/(tod_tumor[time_min_health])))

#     fold_change = tod_healthy/tod_tumor
#     time_max_benefit = toD[np.argmax(fold_change)]
#     time_min_benefit = toD[np.argmin(fold_change)]
#     fold_change_max.append(max(fold_change))
#     fold_change_min.append(min(fold_change))
#     time_max.append(time_max_benefit)
#     time_min.append(time_min_benefit)

# plt.xticks([0,6,12,18,24])
# plt.ylabel('ToD response')
# plt.xlabel('Time of day [hours]')
# # plt.legend(loc='best')
# # plt.savefig(path+'ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()

# # plt.figure(figsize = (12,10))
# # plt.plot(var_vect, max_benefit, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
# # plt.plot(var_vect, min_benefit, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
# # plt.ylabel('Benefit range of ToD response\n healthy tissue vs tumor')
# # plt.xlabel('Hill coefficient') 
# # plt.legend(loc='best')
# # # plt.savefig(path+'Rel_max_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# # plt.show()

# plt.figure(figsize = (12,10))
# plt.plot(var_vect, fold_change_max, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
# plt.plot(var_vect, fold_change_min, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
# plt.ylabel('Fold-change of ToD response\n healthy tissue vs tumor')
# plt.xlabel('Hill coefficient') 
# plt.legend(loc='best')
# plt.savefig(path+'Fold-change_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111, polar=True)
# cmap = plt.get_cmap('YlGn')
# r = np.ones_like(time_max)
# colors = cmap(var_vect)
# scatter = ax.scatter(time_max, r, c=var_vect, s=200, cmap=cmap)
# ax.grid(False)  # Hide the grid
# ax.set_yticklabels([])  # Hide the radial labels
# hour_labels = [f'{i}h' for i in range(0, 24, 3)] 
# ax.set_xticks(np.linspace(0, 2 * np.pi, len(hour_labels), endpoint=False))
# ax.set_xticklabels(hour_labels)
# cbar = plt.colorbar(scatter, ax=ax, orientation='vertical')
# cbar.set_label('Hill coefficient')
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1)
# plt.title('Time of maximal benefit')
# plt.savefig(path+'Time_max_benefit_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()

#%%Change drug stability
#Circadian parameters tumor
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5
d = 0

#Parameters to change
opt = 'change_drug_stability'
rates = np.arange(1, 48, step=1)#[48,24,12,6]  #
kdecay = []
for i in range(len(rates)):
    if rates[i]!=0:
        kdecay.append(1/rates[i])

colors = plt.cm.tab20c(range(len(kdecay)))

max_benefit = []
min_benefit = []

fold_change_max = []
fold_change_min = []

time_max = []
time_min = []

for ii in range(len(kdecay)):
    tod_tumor = []
    
    d = kdecay[ii]

    c2 = circmod_c(modul_strgth, Amp, toD, T, phi)*c_exp(d,toD)
    for j in range(len(c2)):
        tod_tumor.append(growthcurve(1, t, k, kl, LC50, Sm, c2[j], SC50, h)[-1])
    tod_tumor = tod_tumor/tod_tumor[0]
    
    plt.plot(toD, tod_tumor, linewidth=4, color=colors[ii], label=rates[ii]) #####change label
    
    max_benefit.append((max(tod_healthy)-tod_tumor[time_max_health]))
    min_benefit.append((min(tod_healthy)-tod_tumor[time_min_health]))

    # fold_change_max.append((max(tod_healthy)/(tod_tumor[time_max_health])))
    # fold_change_min.append((min(tod_healthy)/(tod_tumor[time_min_health])))

    fold_change = tod_healthy/tod_tumor
    time_max_benefit = toD[np.argmax(fold_change)]
    time_min_benefit = toD[np.argmin(fold_change)]
    fold_change_max.append(max(fold_change))
    fold_change_min.append(min(fold_change))
    time_max.append(time_max_benefit)
    time_min.append(time_min_benefit)
    
plt.xticks([0,6,12,18,24])
plt.ylabel('ToD response')
plt.xlabel('Time of day [hours]')
# plt.legend(loc='best')
# plt.savefig(path+'ToD_healthy_vs_tumor_'+str(opt)+'.svg')
plt.show()


# plt.figure(figsize = (12,10))
# plt.plot(kdecay, max_benefit, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
# plt.plot(kdecay, min_benefit, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
# plt.ylabel('Benefit range of ToD response\n healthy tissue vs tumor')
# plt.xlabel('Decays [1/hours]') 
# plt.legend(loc='best')
# # plt.savefig(path+'Rel_max_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
# plt.show()

plt.figure(figsize = (12,10))
plt.plot(kdecay, fold_change_max, '-o', color='green', markersize=15, linewidth=5, label='Maximum benefit')
plt.plot(kdecay, fold_change_min, '-o', color='tab:red', markersize=15, linewidth=5, label='Minimum benefit')
plt.ylabel('Fold-change of ToD response\n healthy tissue vs tumor')
plt.xlabel('Decays [1/hours]') 
plt.legend(loc='best')
# plt.savefig(path+'Fold-change_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
cmap = plt.get_cmap('YlGn')
r = np.ones_like(time_max)
colors = cmap(kdecay)
scatter = ax.scatter(time_max, r, c=rates, s=200, cmap=cmap)
ax.grid(False)  # Hide the grid
ax.set_yticklabels([])  # Hide the radial labels
hour_labels = [f'{i}h' for i in range(0, 24, 3)] 
ax.set_xticks(np.linspace(0, 2 * np.pi, len(hour_labels), endpoint=False))
ax.set_xticklabels(hour_labels)
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical')
cbar.set_label('Decay rates [hours]')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
plt.title('Time of maximal benefit')
# plt.savefig(path+'Time_max_benefit_ToD_healthy_vs_tumor_'+str(opt)+'.svg')
plt.show()



