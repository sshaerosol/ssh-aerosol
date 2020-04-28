#!/usr/bin/env python
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
pcase_dyn = '../results/cond-evap-inorg/cond-evap-dyn'
pcase_eq = '../results/cond-evap-inorg/cond-evap-eq'
pcase_icut = '../results/cond-evap-inorg/cond-evap-icut'
tag_fig = 'time'
sizebin_ssh = 15
initial_time = 0.0        
final_time = 1000.0             
delta_t = 10.0                 
#####################################################################
font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

nt = int((final_time-initial_time)/delta_t)
conc_nh3_dyn = np.zeros(nt)
conc_nh3_eq = np.zeros(nt)
conc_nh3_icut = np.zeros(nt)
conc_hno3_dyn = np.zeros(nt)
conc_hno3_eq = np.zeros(nt)
conc_hno3_icut = np.zeros(nt)
ssh_time = np.zeros(nt)
for i in range(nt):
   ssh_time[i] = initial_time + delta_t * i 
############################## results from SIREAM ####################################################""
#### 
fig = plt.figure(1,figsize = (15,15))
with open (pcase_dyn+'/gas/NH3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_nh3_dyn[i] = float(values[i])
with open (pcase_eq+'/gas/NH3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_nh3_eq[i] = float(values[i])
with open (pcase_icut+'/gas/NH3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_nh3_icut[i] = float(values[i])
plt.plot(ssh_time, conc_nh3_dyn, label = 'SSH dyn')
plt.plot(ssh_time, conc_nh3_eq, linestyle='-.',label = 'SSH eq')
plt.plot(ssh_time, conc_nh3_icut, linestyle='--',label = 'SSH hyb')

plt.xlabel(r'time (s)')
plt.title( 'NH$_3$ time evolution')
plt.ylabel(r'Concentration ($\mu$g m$^{-3}$)')
plt.legend(loc ='best')		# show legend
fig.savefig('NH3_COND-EVAP_'+tag_fig)

fig.clf()
with open (pcase_dyn+'/gas/HNO3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_hno3_dyn[i] = float(values[i])
with open (pcase_eq+'/gas/HNO3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_hno3_eq[i] = float(values[i])
with open (pcase_icut+'/gas/HNO3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_hno3_icut[i] = float(values[i])
plt.plot(ssh_time, conc_hno3_dyn, label = 'SSH dyn')
plt.plot(ssh_time, conc_hno3_eq, linestyle='-.',label = 'SSH eq')
plt.plot(ssh_time, conc_hno3_icut, linestyle = '--',label = 'SSH hyb')

plt.xlabel(r'time (s)')
plt.title( 'HNO$_3$ time evolution')
plt.ylabel(r'Concentration ($\mu$g m$^{-3}$)')
plt.legend(loc ='best')		# show legend
plt.tight_layout()
fig.savefig('HNO3_COND-EVAP_'+tag_fig)
