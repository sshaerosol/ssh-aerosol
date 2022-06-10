#!/usr/bin/env python3
# -*- coding: iso-8859-1  -*-
import matplotlib
matplotlib.use('Agg')
import argparse
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()

font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

################### input parameters #################################"
filename = ['../results/mcm-wFGL-tag_RO2-0', '../results/mcm-wFGL-tag_RO2-1', '../results/mcm-wFGL-tag_RO2-2', '../results/mcm-wFGL-tag_RO2-3/']
stl = ['b-','g--','k:','r-']
tiny=1E-19
#####################################################################

fig = plt.figure()
for ind,fn in enumerate(filename):
    with open(fn+'/aero/Organic.txt') as f: soa = f.read().splitlines()
    with open(fn+'/gas/BCARY.txt') as f: voc = f.read().splitlines()
    nt = len(soa)
    data = np.zeros(nt)
    voc0,soa0=float(voc[0]),float(soa[0])
    for i in range(nt): 
        data[i] = (float(soa[i])-soa0) /(voc0-float(voc[i])+tiny)*100.
    plt.plot(data,stl[ind],linewidth=1.3,label=fn)
plt.xlabel('time (hour)')
plt.ylabel(r'Concentrations ($\mu$g/m$^3$)')
plt.xlim(0,24)
plt.title('SOA yields (%)')
plt.legend(['0-no RO2','1-generated RO2','2-background RO2','3-generated+background RO2'],loc='best',framealpha=0.5)
fig.savefig('mcm_diffRO2.png')

