#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/toluene-aft/toluene-aft-ideal/aero/'
dir2 = '../results/toluene-aft/toluene-aft-ideal-rdc/aero/'
dir3 = '../results/toluene-aft/toluene-aft-ideal-rdcRm1/aero/'
dir4 = '../results/toluene-aft/toluene-aft-ideal-rdcH2O/aero/'

Nt= 156
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
ide=np.zeros(shape)
uni=np.zeros(shape)
aio=np.zeros(shape)
uniRm1=np.zeros(shape)
uniH2O=np.zeros(shape)


with open (dir1+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  ide[j] = float(string_conc1[j])


with open (dir2+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  uni[j] = float(string_conc1[j])


with open (dir2+'PC7H9O9_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  aio[j] = float(string_conc1[j])


with open (dir3+'Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  uniRm1[j] = float(string_conc1[j])

with open (dir4+'Organic.5.txt') as finit :
        string_conc = finit.read().splitlines()
        for j in range(Nt) :
                  uniH2O[j] = float(string_conc[j])


time=range(5,(Nt+1)*5,5)
time=numpy.array(time)
time=time/60.

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
#plt.ylim(0.,30.)
plt.grid()
plt.plot(time,ide,'k-',label='Ideal')
plt.plot(time,uni,'b--',label='Ideal rdc')
plt.plot(time,aio,'r-.',label='ipso-BPR SOA')
plt.plot(time,uniRm1,'g-',label='Ideal rdc MCM')
plt.plot(time,uniH2O,'m-',label='Ideal rdc H$^2$O')
plt.ylabel("SOA [$\mu g.m^{-3}$]", fontsize=18)
plt.xlabel("time [min]", fontsize=18)
plt.legend(loc=2)
plt.savefig("toluene-aft-rdc.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
