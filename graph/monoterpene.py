#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/mt-ref/'
# '../results/mt-ref/aero/'
dir2 = '../results/mt-rdc/'
# '../results/mt-ref/aero/'

fname = 'concs.txt'
# fname = 'Organic.5.txt

Nt= 1000
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
ref=np.zeros(shape)
rdc=np.zeros(shape)


with open (dir1 + fname) as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
            ref[j] = float(string_conc1[j])

with open (dir2 + fname) as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
            rdc[j] = float(string_conc1[j])
             

dt = 0.1 # in seconds
ratio = 1.0 / dt

time=[x / ratio for x in range(1, (Nt+1)*1, 1)]
time=numpy.array(time)
time=time/60. # conversion to minutes


fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,20.)
plt.grid()
plt.plot(100.0/60., 17.01137, 'ro', label = 'Xavier et al., 2019')
plt.plot(time,ref,'k-',label='Ref')
plt.plot(time,rdc,'b--',label='Rdc')
plt.ylabel("SOA [$\mu g.m^{-3}$]", fontsize=18)
plt.xlabel("time [min]", fontsize=18)
plt.legend(loc=2)
plt.savefig("monoterpene.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
