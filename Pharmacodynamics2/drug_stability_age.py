#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:50:12 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

path ='output/'

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def GRinh(Sm, c, h, SC50): #normalized growth rate inhibition
    return 2**(1-(Sm*c**h)/(SC50**h+c**h)-(1/k)*(kl*c**h)/(LC50**h+c**h))-1

def circmod_c(modul_strgth, Amp, t, T, phi): #circadian modulation complete
    return SC50*(1+modul_strgth*Amp*np.sin(t*(2*np.pi/T)+phi))

def c_exp(kdecay, t): #exponential decay (drug stability)
    return np.exp(-kdecay*t)

def GR(c, inf, EC50, hGR):
    return inf+(1-inf)/(1+(c/EC50)**hGR)

#Simulation parameters
tf = 120
t = np.linspace(0, tf, tf)
    
#Growth parameters
k = np.log(2)/20 #growth rate
doubling_times = np.arange(12, 48, step=1)

#Drug parameters
Sm = 1 #maximal inhibitor effect
SC50 = 1 #half maximum effect
h = 1 #Hill coeff
kl = 0. #killing effect
LC50 = 1 #half killing effect
x0 = 1 

#Circadian parameters
Amp = 1
T = 24
phi = 0 #0 or pi
modul_strgth = 0.5

#Scenario3: Drug age
step = 8
rates = np.arange(0, 48, step=2)#[0,48,24,12,6]
kdecay = []
for i in range(len(rates)):
    if rates[i]!=0:
        kdecay.append(1/rates[i])
    elif rates[i]==0:
        kdecay.append(0)

print(kdecay)

num = 4

#%%Concentration with circadian modulation and decay
c = 1 #EC50
toD = np.arange(0, 25, step=0.1) #time of day
c = c*circmod_c(modul_strgth, Amp, toD, T, phi)*c_exp(kdecay[num],toD)
colors2 = plt.cm.tab20b(np.linspace(0, 1, len(c)))

fig = plt.figure(figsize = (10,10))
time_decay = np.arange(0, 48, step=1)
for i in range(len(kdecay)):
    plt.plot(time_decay, c_exp(kdecay[i], time_decay), label=str(rates[i]))
plt.xticks([0,6,12,18,24,30,36,42])
plt.xlabel('Time [hours]')
plt.ylabel('Concentration')
plt.legend(loc='best')
# plt.savefig(path+'kdecay_vs_time.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
for i in range(len(kdecay)):
    drc = []
    c = circmod_c(modul_strgth, Amp, toD, T, phi)*c_exp(kdecay[i],toD)
    plt.plot(toD, c, label=str(rates[i])+'h')
plt.xlabel('Time')
plt.ylabel('Modulated concentration')
plt.legend(loc='upper left', title='Decay rates')
# plt.savefig(path+'Modulat_concentr_diff_kdecay.svg')
plt.show()

#%%Cell counts 
# fig = plt.figure(figsize = (10,10))
# for j in range(len(c)):
#     plt.plot(t, growthcurve(x0, t, k, kl, LC50, Sm, c[j], SC50, h), color=colors2[j])
#     if j%2==0:
#         plt.plot(t, growthcurve(x0, t, k, kl, LC50, Sm, c[j], SC50, h), color=colors2[j], label=str(int(toD[j]))+'h')
# plt.xlabel('Time')
# plt.ylabel('Normalized cell count')
# plt.legend(loc='upper left', title='Time of day')
# plt.title('Decay rate '+str(rates[num])+'h')
# # plt.savefig(path+'Norm_count_modulat_concentr_diff_kdecay.svg')
# plt.show()

#%%ToD response
mr = []
fig = plt.figure(figsize = (10,10))
for i in range(len(kdecay)):
    drc = []
    c = circmod_c(modul_strgth, Amp, toD, T, phi)*c_exp(kdecay[i],toD)
    for j in range(len(c)):
        ratio = np.log(growthcurve(x0, t, k, kl, LC50, Sm, c[j], SC50, h)[-1])#/growthcurve(x0, t, k, kl, LC50, Sm, c[j], SC50, h)[0]))
        drc.append(ratio)
    mr.append(abs(max(drc/drc[0])-min(drc/drc[0])))
    plt.plot(toD, drc/drc[0], label=str(rates[i])+'h')
plt.xlabel('Time')
plt.ylabel('Time of day response')
plt.legend(loc='upper left', title='Decay rates')
# plt.savefig(path+'ToD_diff_kdecay.svg')
plt.show()

#%%Maximum range of ToD
plt.figure(figsize = (12,10))
plt.plot(kdecay, mr, 'o', markersize=15, linewidth=5)
plt.ylabel('Maximum range relative ToD response')
plt.xlabel('1/decay [1/hours]') 
# plt.savefig(path+'Max_ToD_decays.svg')
plt.show()



