#!/usr/bin/env python3
# -*- coding: iso-8859-1  -*-
import matplotlib
matplotlib.use('Agg')
import math,sys, optparse, glob, os, re, string
from datetime import *
from math import *
from pylab import *
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()

################### input parameters #################################"
pcase1= '../results/condkelv/inorg-dyn/'
pcase2 = '../results/condkelv/inorg-nokelv-dyn/'
pcase3 = '../results/condkelv/inorg-eq/'
pcase4 = '../results/condkelv/inorg-nokelv-eq/'
pcase5= '../results/condkelv/inorg-eq-soapinorg/'
pcase6 = '../results/condkelv/inorg-dyn-soapinorg/'
sizebin_ssh = 50
density = 1.4 * 1e-6 #in ug um-3
iredist = 0 #0 without redistribution; 1 with redistribution
#####################################################################

font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

deltalogd = []
diam = []
sizebin = 50
tmp1 = -1
tmp2 = -1
with open ('cond_ref/number.init') as f1 :
        values = f1.read().splitlines()
for i in range(len(values)) :
        if i <= (sizebin + 2) : continue # sizebound + title + star line
        tmp2 = values[i].split('   ', -1)[0]
        if tmp1 == tmp2 : continue
        else :
                diam.append(tmp2.replace('D','E'))
                tmp1 = tmp2
diam.append(10.)
diam_mean=[]
#### diameter
for j in range(len(diam)-1) :
        diam_mean.append((float(diam[j])*float(diam[j+1]))**0.5)
        deltalogd.append(math.log10(float(diam[j+1]))-math.log10(float(diam[j])))

############################## extract results from simulation
#### 
cases = ['Kelvin Dynamic','No Kelvin Dynamic','Kelv Equilibrium','No Kelvin Equilibrium','Kelvin Equilibrium soap_inorg=1','Kelvin Dynamic soap_inorg=1']
case_lb = ['Kelvin Dynamic','No Kelvin Dynamic','Kelvin Equilibrium','No Kelvin Equilibrium','Kelvin Equilibrium soap_inorg=1','Kelvin Dynamic soap_inorg=1'] 
#### get conc.
num_init_sml = np.zeros((len(cases),sizebin_ssh))
num_out_sml = np.zeros((len(cases),sizebin_ssh))
mass_init_sml = np.zeros((len(cases),sizebin_ssh))
mass_out_sml = np.zeros((len(cases),sizebin_ssh))
cell_diam_init_sml = np.zeros((len(cases),sizebin_ssh))
cell_diam_out_sml = np.zeros((len(cases),sizebin_ssh))

for i in cases :
        for j in range(sizebin_ssh) :
                if (cases.index(i) == 0):
                        pcase = pcase1
                elif (cases.index(i) == 1):
                         pcase = pcase2
                elif (cases.index(i) == 2):
                        pcase = pcase3
                elif (cases.index(i) == 3):
                        pcase = pcase4
                elif (cases.index(i) == 4):
                        pcase = pcase5
                elif (cases.index(i) == 5):
                        pcase = pcase6
                with open (pcase+'/number/NUMBER_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                num_init_sml[cases.index(i)][j] = float(values[0])
                num_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase+'/aero/PSO4_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                mass_init_sml[cases.index(i)][j] = float(values[0])
                mass_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase+'/aero/PNH4_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                mass_init_sml[cases.index(i)][j] += float(values[0])
                mass_out_sml[cases.index(i)][j] += float(values[-1])
                with open (pcase+'/diameter/DIAMETER_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                cell_diam_init_sml[cases.index(i)][j] = float(values[0])
                cell_diam_out_sml[cases.index(i)][j] = float(values[-1])

############################### drawing Particle Number Distribution ############################### 
lists = []
lbs = []
stl = ['-','-','-','-','-']
cols = ['c-+','m-+','b-','g-+','r-+','y-']
for i in num_out_sml : lists.append(i)
for i in case_lb : lbs.append(i)

fig = plt.figure(1,figsize = (15,15))
num = len(deltalogd)
tmp = np.zeros(num)
for j in range(num) :
        tmp[j]= float(num_init_sml[0][j]) / deltalogd[j]* 1E-6 #in cm-3
plt.plot(diam_mean, tmp,'k-+',label = 'Init')

for i in range(len(lists)) :
     if(iredist==1):
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists[i][j]) / deltalogd[j]* 1E-6 #in cm-3
                if (tmp[j] < 0.0001) : tmp[j] = 0.0
        plt.plot(diam_mean, tmp, cols[i],label = lbs[i])
     else:
        tmp = np.zeros(num)
        numtot = 0.0
        for j in range(num) :
                dd = 1.
                if ((j < num-1) and float(cell_diam_out_sml[i][j]) > 0.001):
                   dd = math.log10(float(cell_diam_out_sml[i][j+1]))-math.log10(float(cell_diam_out_sml[i][j]))
                if (j == num-1):
                   dd = math.log10(float(cell_diam_out_sml[i][j]))-math.log10(float(cell_diam_out_sml[i][j-1]))
                if(dd <=0.0):
                  dd = 1.
               # print(cell_diam_out_sml[1][j],cell_diam_out_sml[2][j],cell_diam_out_sml[3][j])
                if(float(cell_diam_out_sml[i][j]) < 0.1):
                    numtot = numtot + float(lists[i][j]) * 1E-6
                tmp[j]= float(lists[i][j]) / dd* 1E-6 #in cm-3
                if (tmp[j] < 0.0001) : 
                        tmp[j] = 0.0
                        cell_diam_out_sml[i][j] = diam_mean[j]
        plt.plot(cell_diam_out_sml[i], tmp, cols[i],label = lbs[i])
#        print (numtot,lbs[i])
#print ('-------------------------------------------------------------')
plt.xlabel(r'd ($\mu$m)')
plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
plt.yscale('log')
plt.title( 'Particle Number Distribution')
plt.legend(loc ='best',fontsize=10)        	# show legend
plt.tight_layout()
fig.savefig('dNdlogd_CONDKELV_INORG')

############################### drawing Particle volume Distribution
lists2 = []
#for i in mass_init_sml : lists2.append(i)
for i in mass_out_sml : lists2.append(i)
lbs = []
cols = ['c-+','m-+','b-','g-+','r-+','y-']
for i in case_lb : lbs.append(i)
plt.clf()
num = len(deltalogd)
tmp = np.zeros(num)
for j in range(num) :
        tmp[j]= float(mass_init_sml[0][j]) / deltalogd[j]/density
        if(tmp[j] < 0.01) : tmp[j] = 0.0
        else:
                tmp[j] = tmp[j] * 1E-6 #in cm-3
plt.plot(diam_mean, tmp,'k-+',label = 'Init')

for i in range(len(lists2)) :
     if(iredist==1):
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists2[i][j])/deltalogd[j]/density
                if(tmp[j] < 0.01) : tmp[j] = 0.0
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

     else:
        tmp = np.zeros(num)
        masstot = 0.0
        for j in range(num) :
                dd = 1.
                if ((j < num-1) and float(cell_diam_out_sml[i][j]) > 0.001):
                   dd = math.log10(float(cell_diam_out_sml[i][j+1]))-math.log10(float(cell_diam_out_sml[i][j]))
                if (j == num-1):
                   dd = math.log10(float(cell_diam_out_sml[i][j]))-math.log10(float(cell_diam_out_sml[i][j-1]))
                if(dd <=0.0):
                  dd = 1.
                tmp[j]= float(lists2[i][j]) / dd /density
                if(float(cell_diam_out_sml[i][j]) < 0.1):
                    masstot = masstot + float(lists2[i][j]) * 1E-6
                if (tmp[j] < 0.01) : 
                        tmp[j] = 0.0
                        cell_diam_out_sml[i][j] = diam_mean[j]
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        plt.plot(cell_diam_out_sml[i], tmp, cols[i],label = lbs[i])
        #print (masstot,lbs[i])

plt.xlabel(r'd ($\mu$m)')
plt.ylabel(r'dV/d log d ($\mu$m$^3$ cm$^{-3}$)')
plt.xscale('log')
plt.title( 'Particle Volume Distribution')
plt.legend(loc ='upper left',fontsize=10)        	# show legend
plt.tight_layout()
fig.savefig('dVdlogd_CONDKELV_INORG')


