#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import matplotlib
matplotlib.use('Agg')
import datetime
from pylab import *
from numpy import *
import numpy as np
from math import pow
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import LogLocator, LogFormatter 
from matplotlib import colors 
from matplotlib import ticker
from PIL import Image

import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
pcase = '../results/cond-ext/'
sizebin_ssh = 50
density = 1.84 * 1e-6 #in ug um-3
N_fraction = 10
composition_bounds = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1.0])

#####################################################################
############################## results from SIREAM ####################################################""
font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

sizebin = 50
diam = []
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
diam.append(10.0)
diam_bounds = np.zeros(sizebin_ssh+1)
#### diameter
for j in range(len(diam)) :
        diam_bounds[j] = diam[j]

num_init = np.zeros((N_fraction,sizebin_ssh))
num_out = np.zeros((N_fraction,sizebin_ssh))
for j in range(sizebin_ssh) :
          for k in range(N_fraction):
                with open (pcase+'/aero/PSO4_' + str(j*N_fraction+(k)+1) + '.txt') as finit :
                     values = finit.read().splitlines()
                num_init[k][j] += float(values[0])
                num_out[k][j] += float(values[-1])
                with open (pcase+'/aero/PMonomer_' + str(j*N_fraction+(k)+1) + '.txt') as finit :
                     values = finit.read().splitlines()
                num_init[k][j] += float(values[0])
                num_out[k][j] += float(values[-1])
print (num_out[0][49], num_init[9][49])

#jetcmap = cm.get_cmap("jet", 500) #generate a jet map with 10 values
#jet_vals = jetcmap(np.arange(500)) #extract those values as an array
#cmap = ListedColormap(jet_vals)
cmap = plt.get_cmap('PuBuGn')

composition_bounds = composition_bounds.transpose()
vmax0 = num_out.max()

figure()
pcolor(diam_bounds,composition_bounds,num_init,cmap=cmap,vmax=vmax0)
ylim(ymin = 0.0, ymax = 1.0)
xticks(diam_bounds)
yticks(composition_bounds)
axes().set_xscale("log")
xlabel(u"diameter [µm]",size=14) #3, fontweight = "bold")
ylabel(u"Fraction of sulfate", size = 14)#, fontweight = "bold")
colorbar()
xlim(xmin = 0.001, xmax = 10)
yticks(composition_bounds)
title(u"Particle mass concentration [$\mu$g m$^{-3}$]",
      fontsize = 14)#, fontweight = "bold")
plt.tight_layout()
savefig('cond_mass_ext_init.png')

figure()
pcolor(diam_bounds,composition_bounds,num_out,cmap=cmap,vmax=vmax0)
ylim(ymin = 0.0, ymax = 1.0)
xticks(diam_bounds)
yticks(composition_bounds)
axes().set_xscale("log")
xlabel(u"diameter [µm]",size=14)#3, fontweight = "bold")
ylabel(u"Fraction of sulfate", size = 14)#, fontweight = "bold")
colorbar()
xlim(xmin = 0.001, xmax = 10)
yticks(composition_bounds)
title(u"Particle mass concentration [$\mu$g m$^{-3}$]",
      fontsize = 14)#, fontweight = "bold")
plt.tight_layout()
savefig('cond_mass_ext_out.png')
