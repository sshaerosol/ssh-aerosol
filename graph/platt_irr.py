#!/usr/bin/env python3
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/aerosol-platt-nh3-isopa/aero/'
dir2 = '../results/aerosol-platt-nh3-isopa-irr/aero/'
dir3 = '../results/aerosol-platt-nh3-isopa-irr-eq/aero/'
#dir4 = '../results/aerosol-platt-nh3-isopa-oligo-inorg/aero/'


Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
org1=np.zeros(shape)
org2=np.zeros(shape)
org3=np.zeros(shape)

irr1=np.zeros(shape)
irr2=np.zeros(shape)
irr3=np.zeros(shape)

with open (dir1+'/Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  org1[j] = org1[j] + float(string_conc1[j])

with open (dir2+'/Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  org2[j] = org2[j] + float(string_conc1[j])
with open (dir2+'/PirrSOA_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  irr2[j] = irr2[j] + float(string_conc1[j])                                    

with open (dir3+'/Organic.5.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  org3[j] = org3[j] + float(string_conc1[j])
with open (dir3+'/PirrSOA_1.txt') as finit :
        string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
                  irr3[j] = irr3[j] + float(string_conc1[j])                                                                                      

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,300.)

plt.plot(time, org1,'k',label='No irreversible SOA')
plt.plot(time, org2,'b',label='With irreversible SOA (Dyn)')
plt.plot(time, irr2,'--b',label='Irreversible SOA (Dyn)')
plt.plot(time, org3,'r',label='With irreversible SOA (Eq)')
plt.plot(time, irr3,'--r',label='Irreversible SOA (Eq)')

plt.ylabel("Isoprene SOA ($\mu g.m^{-3}$)", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=4)
plt.savefig("platt-irreversible.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
