#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir2 = '../results/aerosol-vocox/ref/aero/'
dir3 = '../results/aerosol-vocox/ref-dNO2/aero/'
dir4 = '../results/aerosol-vocox/ref-mNO2/aero/'
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

species_tol={'PBiBmP','PBiBlP'}
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

species_tol={'PC141CO2H','PBCKSOZ','PC132OOH','PC131CO2H','PC133NO3','PC133CO'}
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
plt.plot(time, tol,color='k',label='H$^2$O - Ref')
plt.plot(time, tolmNO2,'m',label='H$^2$O - NO$_2$ x 2')
plt.plot(time, toldNO2,'b',label='H$^2$O - NO$_2$ / 2')
plt.plot(time, tolcompl,color='k', linestyle='dashed',label='Expl - Ref')
plt.plot(time, tolcomplmNO2,'m',linestyle='dashed',label='Expl - NO$_2$ x 2')
plt.plot(time, tolcompldNO2,'b',linestyle='dashed',label='Expl - NO$_2$ / 2')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("particles-vocox-sesqui.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)

