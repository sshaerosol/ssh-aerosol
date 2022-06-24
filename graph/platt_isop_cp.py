#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/aerosol-platt/aero/'
dir2 = '../results/aerosol-platt-nh3-isop/aero/'
dir3 = '../results/aerosol-platt-nh3-isop-cpeq/aero/'
dir4 = '../results/aerosol-platt-nh3-isop-cpdyn/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
inorg=np.zeros(shape)
isop1=np.zeros(shape)
isop2=np.zeros(shape)
isop3=np.zeros(shape)
isop4=np.zeros(shape)

with open (dir1+'/PBiPER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiDER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop1[j] = isop1[j] + float(string_conc1[j])
with open (dir1+'/PBiMT_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop1[j] = isop1[j] + float(string_conc1[j])

with open (dir2+'/PBiPER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiDER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop2[j] = isop2[j] + float(string_conc1[j])
with open (dir2+'/PBiMT_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop2[j] = isop2[j] + float(string_conc1[j])

with open (dir3+'/PBiPER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop3[j] = isop3[j] + float(string_conc1[j])
with open (dir3+'/PBiDER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop3[j] = isop3[j] + float(string_conc1[j])
with open (dir3+'/PBiMT_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop3[j] = isop3[j] + float(string_conc1[j])

with open (dir4+'/PBiPER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop4[j] = isop4[j] + float(string_conc1[j])
with open (dir4+'/PBiDER_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop4[j] = isop4[j] + float(string_conc1[j])
with open (dir4+'/PBiMT_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  isop4[j] = isop4[j] + float(string_conc1[j])

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,80.)

plt.plot(time, isop1,'k-',label='No NH$_3$')
plt.plot(time, isop2,'g-',label='With NH$_3$ - not coupled phases')
plt.plot(time, isop3,'b--',label='With NH$_3$ - coupled phases eq')
plt.plot(time, isop4,'c-.',label='With NH$_3$ - coupled phases dyn')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("platt-isop-cp.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
