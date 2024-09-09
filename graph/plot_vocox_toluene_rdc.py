#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir2 = '../results/aerosol-vocox/ref-rdc/aero/'
dir3 = '../results/aerosol-vocox/ref-rdc-dNO2/aero/'
dir4 = '../results/aerosol-vocox/ref-rdc-mNO2/aero/'
dir5 = '../results/aerosol-vocox/ref-expl/aero/'
dir6 = '../results/aerosol-vocox/ref-expl-dNO2/aero/'
dir7 = '../results/aerosol-vocox/ref-expl-mNO2/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
tol = np.zeros(shape)
toldNO2 = np.zeros(shape)
tolmNO2 = np.zeros(shape)
tolcompl = np.zeros(shape)
tolcompldNO2 = np.zeros(shape)
tolcomplmNO2 = np.zeros(shape)

species_tol= {'PTOL3OH','PAU50DN','PHOM2ONO2','PUD7000','PTOL4OH','PTOL3OH1NO2','PHOM2O','PIRDK3000','PPP4000','PC7H9O9','PMBQN1OH','PMBQN2OH','PMBQN3OH'}
for species in species_tol:
   with open (dir2+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tol[j] = tol[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir3+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  toldNO2[j] = toldNO2[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir4+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolmNO2[j] = tolmNO2[j] + float(string_conc1[j])

species_tol={'PC7H9O9','PP02000','PDK3000','PFURON','PUD4000','PDD3001','PMALAHY','PMFUR','PUD5000','PUD5001','PAR0010','PDK4000','PAR0088','PAR0043','PAU5002','PPU5000','PPU5001','PPU5002','PAR0130','PAR0127','PPK5001','PPP4004','PAR0134','PAR0144','PAR0087','PAR0027','PAR0093','PPK5003','PAR0113','PAR0090','PAR0104','PAR0105','PAR0109','PAR0128','PAR0124','PPH5002','PAR0110','PAR0048','PAR0115','PAR0153','PAR0094','PPH5004','PPP4000','PPH5000','PGH5002','PNTOL','PUD5002','PAU5000','PAU4000','PAD2000','PAR0140OOH','PMBQN1K1OH','PC73K1OHOOH','PDD5002','PMBQN1OH','PTOL3OH','PMBQN2OH','PTOL4OH','PMBQN3OH','PTOL5OH','PTOL3OH1NO2','PTOL2OHOOH','PTOL4OH1NO2','PTOL3OHOOH','PBTOL3OHOOH','PAD5000','PDK5000','PAD4000','PAA2000','PBTOL4OHOOH','PAD4001','PAK5000','PBTOL5OHOOH','PAA4000','PUD6000','PUD7000','PAR0132','PBZALDOH','PMe6Cy1U3K','PFURR3','PFURR5','PAU6000','PDD3000','PDK4001','PAU7000','PUU7000','PAK3000','PFURR6','PFUROH','PA02000','PHOM1O','PHOM1ONO2','PHOM1OOH','PHOM2O','PHOM2ONO2','PHOM2OOH','PFURR6OHOOH','PAU50DN','PAR0138','PHD2000','PED5000','PED5000OOH','PED5002OOH','PED4001','PIRDK3000'}
for species in species_tol:
   with open (dir5+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcompl[j] = tolcompl[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir6+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcompldNO2[j] = tolcompldNO2[j] + float(string_conc1[j])

for species in species_tol:
   with open (dir7+'/'+species+'_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  tolcomplmNO2[j] = tolcomplmNO2[j] + float(string_conc1[j])


time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,1.)

#plt.plot(time, totsoa,'k-',label='total')
plt.plot(time, tol,color='k',label='Rdc - Ref')
plt.plot(time, tolmNO2,'m',label='Rdc - NO$_2$ x 1.5')
plt.plot(time, toldNO2,'b',label='Rdc - NO$_2$ / 1.5')
plt.plot(time, tolcompl,color='k', linestyle='dashed',label='Expl. - Ref')
plt.plot(time, tolcomplmNO2,'m',linestyle='dashed',label='Expl. - NO$_2$ x 1.5')
plt.plot(time, tolcompldNO2,'b',linestyle='dashed',label='Expl. - NO$_2$ / 1.5')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("particles-vocox-tol-rdc.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

