#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/naphthalene-aft/naphthalene-aft-ideal/aero/'
dir2 = '../results/naphthalene-aft/naphthalene-aft-ideal-no2/aero/'
dir3 = '../results/naphthalene-aft/naphthalene-aft-aiomfac/aero/'
dir4 = '../results/naphthalene-aft/naphthalene-aft-aiomfac-no2/aero/'
dir5 = '../results/naphthalene-aft/naphthalene-aft-aiomfac_visc/aero/'
dir6 = '../results/naphthalene-aft/naphthalene-aft-aiomfac_visc-no2/aero/'

Nt= 156
#####################################################################


font = {'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=12)

shape = [Nt]
ide=np.zeros(shape)
ide_no2=np.zeros(shape)
aio=np.zeros(shape)
aio_no2=np.zeros(shape)
aio_visc=np.zeros(shape)
aio_visc_no2=np.zeros(shape)

with open (dir1+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  ide[j] = float(string_conc1[j])
                  
with open (dir2+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  ide_no2[j] = float(string_conc1[j])


with open (dir3+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  aio[j] = float(string_conc1[j])
                  
with open (dir4+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  aio_no2[j] = float(string_conc1[j])
                  
                  
with open (dir5+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  aio_visc[j] = float(string_conc1[j])

with open (dir6+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  aio_visc_no2[j] = float(string_conc1[j])
                                   

time=range(5,(Nt+1)*5,5)
time=numpy.array(time)
time=time/60.

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,20.)
plt.grid()
plt.plot(time,ide,'k-',label='Ideal')
plt.plot(time,ide_no2,'k--',label='Ideal + 200ppb NO2')
plt.plot(time,aio,'b-',label='AIOMFAC')
plt.plot(time,aio_no2,'b--',label='AIOMFAC + 200ppb NO2')
plt.plot(time,aio_visc,'r-',label='AIOMFAC-visc')
plt.plot(time,aio_visc_no2,'r--',label='AIOMFAC-visc + 200ppb NO2')
plt.ylabel("SOA [$\mu g.m^{-3}$]", fontsize=18)
plt.xlabel("time [min]", fontsize=18)
plt.legend(loc=2)
plt.savefig("naphthalene-aft.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
