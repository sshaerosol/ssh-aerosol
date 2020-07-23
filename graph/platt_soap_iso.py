#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir1 = '../results/aerosol-platt-nh3-isopu_soap/aero/'
dir2 = '../results/aerosol-platt-nh3-isopu_iso/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
nit1=np.zeros(shape)
nit2=np.zeros(shape)
wat1=np.zeros(shape)
wat2=np.zeros(shape)
amm1=np.zeros(shape)
amm2=np.zeros(shape)

with open (dir1+'/PNO3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   nit1[j] = nit1[j] + float(string_conc1[j])
with open (dir2+'/PNO3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   nit2[j] = nit2[j] + float(string_conc1[j])

with open (dir1+'/PNH4_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   amm1[j] = amm1[j] + float(string_conc1[j])
with open (dir2+'/PNH4_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   amm2[j] = amm2[j] + float(string_conc1[j])

with open (dir1+'/PH2O_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   wat1[j] = wat1[j] + float(string_conc1[j])
with open (dir2+'/PH2O_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   wat2[j] = wat2[j] + float(string_conc1[j])

with open (dir1+'/PNO3_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   nit1[j] = nit1[j] + float(string_conc1[j])
with open (dir2+'/PNO3_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   nit2[j] = nit2[j] + float(string_conc1[j])

with open (dir1+'/PNH4_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   amm1[j] = amm1[j] + float(string_conc1[j])
with open (dir2+'/PNH4_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   amm2[j] = amm2[j] + float(string_conc1[j])

with open (dir1+'/PH2O_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   wat1[j] = wat1[j] + float(string_conc1[j])
with open (dir2+'/PH2O_2.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   wat2[j] = wat2[j] + float(string_conc1[j])

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,30.)

plt.plot(time, nit1,'k-',label='NO3 SOAP')
plt.plot(time, nit2,'g-',label='NO3 ISORROPIA')
plt.plot(time, amm1,'k--',label='NH4 SOAP')
plt.plot(time, amm2,'g--',label='NH4 ISORROPIA')
plt.plot(time, wat1,'k-.',label='WATER SOAP')
plt.plot(time, wat2,'g-.',label='WATER ISORROPIA')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc="lower right")
plt.savefig("platt-soap-iso_nit.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
