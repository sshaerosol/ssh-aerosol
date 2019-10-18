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
pcase = '../results/coag/'
sizebin_ssh = 50
density = 1.84 * 1e-6 #in ug um-3

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

############################## results from SIREAM ####################################################""
sizebin = 50
tmp1 = -1
tmp2 = -1
with open ('coag_ref/number.init') as f1 :
	values = f1.read().splitlines()
with open ('coag_ref/number.out') as f2 :
	values2 = f2.read().splitlines()
with open ('coag_ref/comp.init') as ff1 :
	values3 = ff1.read().splitlines()
with open ('coag_ref/comp.out') as ff2 :
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

############################## extract results from Zhang ###########################################
with open ('coag_ref/coag_ni_exact.txt') as fexai :
	for j in range(60):
	    exti = fexai.readline()
            str_exti = exti.split()
	    exact_num_out.append(str_exti[1])
    	    exa_diam_out.append(float(str_exti[0]))
with open ('coag_ref/coag_vi_exact.txt') as fexa2 :
	for j in range(60):
	    exti = fexa2.readline()
            str_exti = exti.split()
	    exact_mass_out.append(float(str_exti[1]))


#### diameter
for j in range(len(exa_diam_out)) :
	if (j == 0) : dt_exa_diam_out.append(math.log10(float(exa_diam_out[j])))
	if (j != 0) : dt_exa_diam_out.append(math.log10(float(exa_diam_out[j]))-math.log10(float(exa_diam_out[j-1])))
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
lists = [num_org_init, num_org_out]
lbs = ['SIREAM_init', 'SIREAM_out']
stl = ['-','-','-','-','-','-','-']
cols = ['*','*','-','-','.','-','.','-']
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
		if (tmp[j] < 0.01) : tmp[j] = 0.0
		else:
			tmp[j] = tmp[j] * 1E-6 #in cm-3
	plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

num = len(exa_diam_out)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(exact_num_out[j])
	if(tmp[j] < 0.01) : tmp[j] = 0.0
#plt.plot(exa_diam_out, tmp,'-',label = 'Zhang_out')
import numpy.ma as ma
tmp = np.array(tmp)
masked = ma.masked_where(tmp == 0.0, tmp)
plt.plot(exa_diam_out, masked,'-',label = 'Zhang_out')


plt.xlabel(r'd($\mu$m)')
plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
plt.yscale('log')
plt.title( 'Particle Number Distribution - case COAG')
plt.legend(loc ='best')		# show legend
fig.savefig('dNdlogd_COAG')

############################### drawing Particle volume Distribution
lists2 = [mass_org_init, mass_org_out]
for i in mass_init_sml : lists2.append(i)
for i in mass_out_sml : lists2.append(i)
lbs = ['SIREAM_init', 'SIREAM_out']
#cols = ['b-','b-','y*','y*','r.','r-','g-','g-']
cols = ['*','*','-','-','-','-','-','-']
for i in case_lb : lbs.append(i + '_init')
for i in case_lb : lbs.append(i + '_out')
# lists = [org_init, org_out, num_init_sml[], num_init_sml[]]
fig = plt.figure(2,figsize = (15,15))
num = len(deltalogd)
for i in range(len(lists2)) :
	tmp = np.zeros(num)
        for j in range(num) :
		tmp[j]= float(lists2[i][j])/deltalogd[j]/density
		if(tmp[j] < 0.01) : tmp[j] = 0.0
		else:
			tmp[j] = tmp[j] * 1E-6 #in cm-3
	plt.plot(diam_mean, tmp, cols[i],label = lbs[i])

num = len(exa_diam_out)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(exact_mass_out[j]) 
plt.plot(exa_diam_out, exact_mass_out,'-',label = 'Zhang_out')

plt.xlabel(r'd($\mu$m)')
#plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
#plt.yscale('log')
plt.title( 'Particle Volume Distribution - case COAG')
plt.legend(loc ='best')		# show legend
fig.savefig('dVdlogd_COAG')


