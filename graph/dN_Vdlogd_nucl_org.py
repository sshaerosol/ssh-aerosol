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
pcasea = '../results/nucl_coupl/'
pcaseb = '../results/nucl_org/'
tag_fig = 'nucl_org'
sizebin_ssh = 50
density = 1.56 * 1e-6 #in ug um-3
iredist = 1 #0 without redistribution; 1 with redistribution
#iredist = 1 #0 without redistribution
#####################################################################
num_org_init = []
num_org_out = []
mass_org_init = []
mass_org_out = []
vol_org_init = []
vol_org_out = []
diam = []
deltalogd = []

exact_num_init = []
exact_mass_init = []
exact_num_out = []
exact_mass_out = []
exa_diam_init = []
dt_exa_diam_init = []
exa_diam_out = []
dt_exa_diam_out = []

font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)
############################## results from SIREAM ####################################################""
sizebin = 50
tmp1 = -1
tmp2 = -1
with open ('cond_ref/number.init') as f1 :
        values = f1.read().splitlines()
with open ('cond_ref/number.out') as f2 :
        values2 = f2.read().splitlines()
with open ('cond_ref/comp.init') as ff1 :
        values3 = ff1.read().splitlines()
with open ('cond_ref/comp.out') as ff2 :
        values4 = ff2.read().splitlines()
for i in range(len(values)) :
        if i <= (sizebin + 2) : continue # sizebound + title + star line
        tmp2 = values[i].split('   ', -1)[0]
        if tmp1 == tmp2 : continue
        else :
                num_org_init.append('%.16E'%(float(values[i].split('   ', -1)[1].replace('D','E'))*1e6)) #in # m-3
                num_org_out.append('%.16E'%(float(values2[i].split('   ', -1)[1].replace('D','E'))*1e6)) #in # m-3
                mass_org_init.append('%.16E'%(float(values3[i].split('   ', -1)[2].replace('D','E')))) # in ug m-3
                mass_org_out.append('%.16E'%(float(values4[i].split('   ', -1)[2].replace('D','E')))) # in ug m-3
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
#cases = ['SSH-coupled','SSH-splitted','SSH-monomer','SSH-monomer-eq']
cases = ['SSH','SSH-monomer']#,'SSH-monomer-dyn']
pcase = [pcasea,pcaseb]#,pcased]
#### get conc.
num_init_sml = np.zeros((len(cases),sizebin_ssh))
num_out_sml = np.zeros((len(cases),sizebin_ssh))
mass_init_sml = np.zeros((len(cases),sizebin_ssh))
mass_out_sml = np.zeros((len(cases),sizebin_ssh))
cell_diam_init_sml = np.zeros(sizebin_ssh)
cell_diam_out_sml = np.zeros(sizebin_ssh)

for i in cases :
        for j in range(sizebin_ssh) :
                with open (pcase[cases.index(i)] +'/number/NUMBER_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                num_init_sml[cases.index(i)][j] = float(values[0])
                num_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase[cases.index(i)] +'/aero/PSO4_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                mass_init_sml[cases.index(i)][j] = float(values[0])
                mass_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase[cases.index(i)] +'/aero/PNH4_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                mass_init_sml[cases.index(i)][j] = float(values[0])
                mass_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase[cases.index(i)] +'/diameter/DIAMETER_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                cell_diam_init_sml[j] = float(values[0])
                cell_diam_out_sml[j] = float(values[-1])

############################### drawing Particle Number Distribution
lists = [num_org_init, num_org_out]
lbs = ['SIREAM_init', 'SIREAM_out']
stl = ['-','-','-','-','-','-','-']
cols = ['-','-','-','-o','-x','-','-','-','-','-']
#for i in num_init_sml : lists.append(i)
lists.append(num_init_sml[cases.index('SSH')])
for i in num_out_sml : lists.append(i)
#for i in cases : lbs.append(i + '_init')
lbs.append('SSH_init')
for i in cases : lbs.append(i + '_out')

fig = plt.figure(1,figsize = (15,15))
num = len(deltalogd)
for i in range(len(lists)) :
     if(i < 3) or (iredist==1):
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists[i][j]) / deltalogd[j]* 1E-6 #in cm-3
                if (tmp[j] < 0.00001) : tmp[j] = 0.0
        if(i>=2):
            tmp = np.array(tmp)
            masked = ma.masked_where(tmp == 0.0, tmp)
            plt.plot(diam_mean, masked, cols[i],label = lbs[i])
     else:
        tmp = np.zeros(num)
        for j in range(num) :
                dd = 1.
                if ((j < num-1) and float(cell_diam_out_sml[j]) > 0.001):
                   dd = math.log10(float(cell_diam_out_sml[j+1]))-math.log10(float(cell_diam_out_sml[j]))
                if (j == num-1):
                   dd = math.log10(float(cell_diam_out_sml[j]))-math.log10(float(cell_diam_out_sml[j-1]))
                if(dd <=0.0):
                  dd = 1.
                tmp[j]= float(lists[i][j]) / dd* 1E-6 #in cm-3
                if (tmp[j] < 0.00001) : tmp[j] = 0.0
        tmp = np.array(tmp)
        masked = ma.masked_where(tmp == 0.0, tmp)
        plt.plot(cell_diam_out_sml, masked, cols[i],label = lbs[i])


plt.xlabel(r'd($\mu$m)')
plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
plt.yscale('log')
plt.title( 'Particle Number Distribution') # - case NUCL')
plt.legend(loc ='best')                # show legend
plt.tight_layout()
fig.savefig('dNdlogd_'+tag_fig)

############################### drawing Particle volume Distribution
lists2 = [mass_org_init, mass_org_out]
#for i in mass_init_sml : lists2.append(i)
lists2.append(mass_init_sml[cases.index('SSH')])
for i in mass_out_sml : lists2.append(i)
#lbs = ['SIREAM_init', 'SIREAM_out']
#for i in cases : lbs.append(i + '_init')
#for i in cases : lbs.append(i + '_out')
# lists = [org_init, org_out, num_init_sml[], num_init_sml[]]
#fig = plt.figure(2,figsize = (15,15))
plt.clf()
num = len(deltalogd)
for i in range(len(lists2)) :
     if(i < 3) or (iredist==1):
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists2[i][j])/deltalogd[j]/density
                if(tmp[j] < 0.01) : tmp[j] = 0.0
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        if(i>=2):
           plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

     else:
        tmp = np.zeros(num)
        for j in range(num) :
                dd = 1.
                if ((j < num-1) and float(cell_diam_out_sml[j]) > 0.001):
                   dd = math.log10(float(cell_diam_out_sml[j+1]))-math.log10(float(cell_diam_out_sml[j]))
                if (j == num-1):
                   dd = math.log10(float(cell_diam_out_sml[j]))-math.log10(float(cell_diam_out_sml[j-1]))
                if(dd <=0.0):
                  dd = 1.
                tmp[j]= float(lists2[i][j]) / dd /density
                if (tmp[j] < 0.01) : tmp[j] = 0.0
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        plt.plot(cell_diam_out_sml, tmp, cols[i],label = lbs[i])

plt.xlabel(r'd($\mu$m)')
plt.ylabel(r'dV/d log d (cm$^{-3}$)')
plt.xscale('log')
#plt.yscale('log')
plt.title( 'Particle Volume Distribution')# - case NUCL')
plt.legend(loc ='best')                # show legend
plt.tight_layout()
fig.savefig('dVdlogd_'+tag_fig)


