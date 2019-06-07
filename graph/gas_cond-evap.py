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
from atmopy import *
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()

################### input parameters #################################"
pcase_dyn = '../results/cond-evap-inorg/cond-evap-dyn'
pcase_eq = '../results/cond-evap-inorg/cond-evap-eq'
tag_fig = 'time'
sizebin_ssh = 15
initial_time = 0.0        
final_time = 1000.0             
delta_t = 10.0                 
#####################################################################

nt = int((final_time-initial_time)/delta_t)
conc_nh3_dyn = np.zeros(nt)
conc_nh3_eq = np.zeros(nt)
conc_hno3_dyn = np.zeros(nt)
conc_hno3_eq = np.zeros(nt)
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
plt.plot(ssh_time, conc_nh3_dyn, label = 'SSH dyn')
plt.plot(ssh_time, conc_nh3_eq, label = 'SSH eq')

plt.xlabel(r'time (s)')
plt.title( 'NH3 time evolution')
plt.legend(loc ='best')		# show legend
fig.savefig('NH3_COND-EVAP_'+tag_fig)


fig = plt.figure(2,figsize = (15,15))
with open (pcase_dyn+'/gas/HNO3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_hno3_dyn[i] = float(values[i])
with open (pcase_eq+'/gas/HNO3.txt') as finit :
	values = finit.read().splitlines()
	for i in range(nt):
           conc_hno3_eq[i] = float(values[i])
plt.plot(ssh_time, conc_hno3_dyn, label = 'SSH dyn')
plt.plot(ssh_time, conc_hno3_eq, label = 'SSH eq')

plt.xlabel(r'time (s)')
plt.title( 'HNO3 time evolution')
plt.legend(loc ='best')		# show legend
fig.savefig('HNO3_COND-EVAP_'+tag_fig)
