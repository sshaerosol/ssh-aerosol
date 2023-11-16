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
#warnings.filterwarnings("ignore",category =RuntimeWarning)

################### input parameters #################################"
pcase = '../results/nucl_org_only_d25/'
sizebin_ssh = 25
density = 1.4 * 1e-6 #in ug um-3

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
####################################################################""
sizebin = 25
diam = [1.00000000e-09,1.44543977e-09,2.08929613e-09,3.01995172e-09,
 4.36515832e-09, 6.30957344e-09, 9.12010839e-09, 1.31825674e-08,
 1.90546072e-08, 2.75422870e-08, 3.98107171e-08, 5.75439937e-08,
 8.31763771e-08, 1.20226443e-07, 1.73780083e-07, 2.51188643e-07,
 3.63078055e-07, 5.24807460e-07, 7.58577575e-07, 1.09647820e-06,
 1.58489319e-06, 2.29086765e-06, 3.31131121e-06, 4.78630092e-06,
 6.91830971e-06, 1.00000000e-05]

diam_mean=[]
#### diameter
for j in range(len(diam)-1) :
        diam_mean.append((float(diam[j])*float(diam[j+1]))**0.5)
        deltalogd.append(math.log10(float(diam[j+1]))-math.log10(float(diam[j])))

############################## extract results from simulation
#### 
cases = ['SSH'] #, 'num1res0']
case_lb = ['SSH'] #, 'num1']
#### get conc.
num_init_sml = np.zeros((len(cases),sizebin_ssh))
num_out_sml = np.zeros((len(cases),sizebin_ssh))
mass_init_sml = np.zeros((len(cases),sizebin_ssh))
mass_out_sml = np.zeros((len(cases),sizebin_ssh))
for i in cases :
        for j in range(sizebin_ssh) :
                with open (pcase+'/number/NUMBER_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                num_init_sml[cases.index(i)][j] = float(values[0])
                num_out_sml[cases.index(i)][j] = float(values[-1])
                with open (pcase+'/aero/PSO4_' + str(j+1) + '.txt') as finit :
                        values = finit.read().splitlines()
                mass_init_sml[cases.index(i)][j] = float(values[0])
                mass_out_sml[cases.index(i)][j] = float(values[-1])

############################### drawing Particle Number Distribution
lists = []
lbs = []
stl = ['-','-','-','-','-']
cols = ['-','-','.','-','.','-']
for i in num_init_sml : lists.append(i)
for i in num_out_sml : lists.append(i)
for i in case_lb : lbs.append(i + '_init')
for i in case_lb : lbs.append(i + '_out')

fig = plt.figure(1,figsize = (15,15))
num = len(deltalogd)
for i in range(len(lists)) :
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists[i][j]) / deltalogd[j]
                if (tmp[j] < 0.1) : tmp[j] = 0.0
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

plt.xlabel(r'd ($\mu$m)')
plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
plt.yscale('log')
plt.title( 'Particle Number Distribution')
plt.legend(loc ='best')        	# show legend
plt.tight_layout()
fig.savefig('dNdlogd_NUCL-ORG_ONLY-d25')

############################### drawing Particle volume Distribution
lists2 = []
for i in mass_init_sml : lists2.append(i)
for i in mass_out_sml : lists2.append(i)
lbs = []
#cols = ['b-','b-','y*','y*','r.','r-','g-','g-']
cols = ['-','-','-','-','-','-']
for i in case_lb : lbs.append(i + '_init')
for i in case_lb : lbs.append(i + '_out')
# lists = [org_init, org_out, num_init_sml[], num_init_sml[]]
#fig = plt.figure(2,figsize = (15,15))
plt.clf()
num = len(deltalogd)
for i in range(len(lists2)) :
        tmp = np.zeros(num)
        for j in range(num) :
                tmp[j]= float(lists2[i][j])/deltalogd[j]/density
                if(tmp[j] < 0.01) : tmp[j] = 0.0
                else:
                        tmp[j] = tmp[j] * 1E-6 #in cm-3
        plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

plt.xlabel(r'd ($\mu$m)')
plt.ylabel(r'dV/d log d ($\mu$m$^3$ cm$^{-3}$)')
plt.xscale('log')
#plt.yscale('log')
plt.title( 'Particle Volume Distribution')
plt.legend(loc ='best')        	# show legend
plt.tight_layout()
fig.savefig('dVdlogd_NUCL-ORG_ONLY-d25')


