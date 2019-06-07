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
pcase_kelv = '../results/kelvin/w_kelv'
pcase_nokelv = '../results/kelvin/no_kelv/'
density = 1.3
iredist = 1 #0 without redistribution; 1 with redistribution
#####################################################################
pcase = [pcase_kelv,pcase_nokelv]
sizebin_ssh = 50
num_org_init = []
num_org_out = []
mass_org_init = []
mass_org_out = []
vol_org_init = []
vol_org_out = []
diam = []
deltalogd = []
deltalogd_out = []
deltalogd_out_1b = []

exact_num_init = []
exact_mass_init = []
exact_num_out = []
exact_mass_out = []
exa_diam_init = []
dt_exa_diam_init = []
exa_diam_out = []
dt_exa_diam_out = []


############################## results from SIREAM COND to read diameters ####################################################""
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

############################## extract results from Zhang ###########################################

############################## results from Devilliers ####################################################""
######################## 500 sections Devilliers #########################

pts        = 'devilliers/points.txt'
dist_init_nb = 'devilliers/dist_init_nb.txt'
dist_init_vol = 'devilliers/dist_init_vol.txt'

nb_full = 'devilliers/graphe_NB_full_moving_sim2.txt'
vol_full = 'devilliers/graphe_VOL_full_moving_sim2.txt'
diam_full = 'devilliers/graphe_diam_full_moving_sim2.txt'


contents_13 = open(pts)
contents_14 = open(dist_init_nb)
contents_15 = open(dist_init_vol)

contents_16 = open(nb_full)
contents_17 = open(vol_full)
contents_18 = open(diam_full)

string_13 = contents_13.readline()
string_14 = contents_14.readline()
string_15 = contents_15.readline()
string_16 = contents_16.readline()
string_17 = contents_17.readline()
string_18 = contents_18.readline()

points   = string_13.split()
dist_nb  = string_14.split()
dist_vol = string_15.split()

full_nb = string_16.split()
full_vol = string_17.split()
full_diam = string_18.split()

num_point = range(600)
num_point_full = range(500)

for j in num_point:
    dist_nb[j]  = float(dist_nb[j])
    dist_vol[j] = float(dist_vol[j])
    points[j]   = float(points[j])

for k in num_point_full:
    full_nb[k] = float(full_nb[k])
    full_vol[k] = float(full_vol[k])
    full_diam[k] = float(full_diam[k])


######################## 48 sections Devilliers Hemen #########################
#/profils_cerea/devillim/General_21_11_11/RESULTS_simulation_2/graphe/48_sect/

nb_sect    = 'devilliers/nb_sections_48sect.txt'
diam_init      = 'devilliers/graphe_diam_init_48sect.txt'
nb_final_e_h  = 'devilliers/graphe_NB_fin_e_h_48sect_sim2.txt'
vol_final_e_h  = 'devilliers/graphe_VOL_fin_e_h_48sect_sim2.txt'

contents_0 = open(nb_sect)
contents_1 = open(diam_init)
contents_4 = open(nb_final_e_h)
contents_9 = open(vol_final_e_h)
ns       = contents_0.readline()
string_1 = contents_1.readline()
string_4 = contents_4.readline()
string_9 = contents_9.readline()
section = arange(int(ns))
dinit      = string_1.split()
nfin_e_h   = string_4.split()
vfin_e_h   = string_9.split()

for i in section:
    dinit[i] = float(dinit[i])
    nfin_e_h[i]  = float(nfin_e_h[i])



############################## extract results from simulation
#### 
cases = ['SSH','SSH_noKelv'] 
#### get conc. Kelvin
num_init_sml = np.zeros((len(cases),sizebin_ssh))
num_out_sml = np.zeros((len(cases),sizebin_ssh))
mass_init_sml = np.zeros((len(cases),sizebin_ssh))
mass_out_sml = np.zeros((len(cases),sizebin_ssh))
cell_diam_init_sml = np.zeros((len(cases),sizebin_ssh))
cell_diam_out_sml = np.zeros((len(cases),sizebin_ssh))
deltalogd_out_sml = np.zeros((len(cases),sizebin_ssh))
deltalogd_init_sml = np.zeros((len(cases),sizebin_ssh))

for i in cases :
	for j in range(sizebin_ssh) :
		with open (pcase[cases.index(i)]+'/number/NUMBER_' + str(j+1) + '.txt') as finit :
			values = finit.read().splitlines()
		num_init_sml[cases.index(i)][j] = float(values[0])
		num_out_sml[cases.index(i)][j] = float(values[-1])
		with open (pcase[cases.index(i)]+'/aero/PPOAmP_' + str(j+1) + '.txt') as finit :
			values = finit.read().splitlines()
		mass_init_sml[cases.index(i)][j] = float(values[0])
		mass_out_sml[cases.index(i)][j] = float(values[-1])
		with open (pcase[cases.index(i)]+'/diameter/DIAMETER_' + str(j+1) + '.txt') as finit :
			values = finit.read().splitlines()
		cell_diam_init_sml[cases.index(i)][j] = float(values[0])
		cell_diam_out_sml[cases.index(i)][j] = float(values[-1])

	for j in range(sizebin_ssh-1):
	    deltalogd_init_sml[cases.index(i)][j] = math.log10(float(cell_diam_init_sml[cases.index(i)][j+1]))-math.log10(float(cell_diam_init_sml[cases.index(i)][j]))
	j = sizebin_ssh-1
	deltalogd_init_sml[cases.index(i)][j] = math.log10(float(cell_diam_init_sml[cases.index(i)][j]))-math.log10(float(cell_diam_init_sml[cases.index(i)][j-1]))

	for j in range(sizebin_ssh-1):
	    deltalogd_out_sml[cases.index(i)][j] = abs(math.log10(float(cell_diam_out_sml[cases.index(i)][j+1]))-math.log10(float(cell_diam_out_sml[cases.index(i)][j])))
	j = sizebin_ssh-1
	deltalogd_out_sml[cases.index(i)][j] =  abs(math.log10(float(cell_diam_out_sml[cases.index(i)][j]))-math.log10(float(cell_diam_out_sml[cases.index(i)][j-1])))

############################### drawing Particle Number Distribution#################
lists = []
lists_diam = []
lists_logdiam = []
lbs=[]
stl = ['-','-','-','-']
cols = ['-','-','-','-']
for i in num_init_sml : lists.append(i)
for i in num_out_sml : lists.append(i)
for i in cell_diam_init_sml : lists_diam.append(i)
for i in cell_diam_out_sml : lists_diam.append(i)
for i in deltalogd_init_sml : lists_logdiam.append(i)
for i in deltalogd_out_sml : lists_logdiam.append(i)
for i in cases : lbs.append(i + '_init')
for i in cases : lbs.append(i + '_out')

fig = plt.figure(1,figsize = (15,15))
num = len(deltalogd)
for i in range(len(lists)) :
   if(i!=1):
     if(i < 1) or (iredist==1):
	tmp = np.zeros(num)
	for j in range(num) :
		tmp[j]= float(lists[i][j]) / deltalogd[j]
		if (tmp[j] < 0.01) : tmp[j] = 0.0
		else:
			tmp[j] = tmp[j] * 1E-6 #in cm-3
	plt.plot(diam_mean, tmp, cols[i],label = lbs[i])
     else:
	tmp = np.zeros(num)
	diam_tmp = np.zeros(num)
	for j in range(num) :
                diam_tmp[j]  =float(lists_diam[i][j])
                if(diam_tmp[j] < 0.001):
			diam_tmp[j]  = 0.001
                        tmp[j] = 0.0
		tmp[j]= float(lists[i][j]) / float(lists_logdiam[i][j])
		if (tmp[j] < 0.01) : tmp[j] = 0.0
		else:
			tmp[j] = tmp[j] * 1E-6 #in cm-3
	plt.plot(diam_tmp, tmp, cols[i],label = lbs[i])

num = len(dist_nb)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(dist_nb[j])# / dt_exa_diam_out[j]
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(points, tmp,'-',label = 'Devilliers_init')

num = len(full_nb)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(full_nb[j])# / dt_exa_diam_out[j]
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(full_diam, tmp,'-',label = 'Devilliers_out')

num = len(nfin_e_h)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(nfin_e_h[j])# / dt_exa_diam_out[j]
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(dinit, tmp,'*',label = 'Devilliers_48bin_out')

plt.xlabel(r'd($\mu$m)')
plt.ylabel(r'dN/d log d (cm$^{-3}$)')
plt.xscale('log')
plt.yscale('log')
plt.title( 'Particle Number Distribution - case KELVIN')
plt.legend(loc ='best')		# show legend
fig.savefig('dNdlogd_KELVIN')

############################### drawing Particle volume Distribution
lists2 = []
for i in mass_init_sml : lists2.append(i)
for i in mass_out_sml : lists2.append(i)
fig = plt.figure(2,figsize = (15,15))
num = len(deltalogd)
for i in range(len(lists)) :
   if(i!=1):
     if(i < 1) or (iredist==1):
	tmp = np.zeros(num)
	for j in range(num) :
		tmp[j]= float(lists2[i][j]) / deltalogd[j]/density  
		if (tmp[j] < 0.01) : tmp[j] = 0.0
	plt.plot(diam_mean, tmp, cols[i],label = lbs[i])
     else:
	tmp = np.zeros(num)
	diam_tmp = np.zeros(num)
	for j in range(num) :
                diam_tmp[j]  =float(lists_diam[i][j])
                if(diam_tmp[j] < 0.001):
			diam_tmp[j]  = 0.001
                        tmp[j] = 0.0
		tmp[j]= float(lists2[i][j])  / float(lists_logdiam[i][j])/density  
		if (tmp[j] < 0.01) : tmp[j] = 0.0
	plt.plot(diam_tmp, tmp, cols[i],label = lbs[i])

num = len(dist_nb)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j] = float(dist_vol[j])
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(points, tmp,'-',label = 'Devilliers_init')

num = len(full_nb)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(full_vol[j])# / dt_exa_diam_out[j]
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(full_diam, tmp,'-',label = 'Devilliers_out')

num = len(vfin_e_h)
tmp = np.zeros(num)
for j in range(num) :
	tmp[j]= float(vfin_e_h[j])# / dt_exa_diam_out[j]
	if(tmp[j] < 0.01) : tmp[j] = 0.0
plt.plot(dinit, tmp,'*',label = 'Devilliers_48bin_out')

plt.xlabel(r'd($\mu$m)')
plt.xscale('log')
plt.title( 'Particle Volume Distribution - case KELVIN')
plt.legend(loc ='best')		# show legend
fig.savefig('dVdlogd_KELVIN')


