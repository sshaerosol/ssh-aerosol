#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy, pylab

# Turn interactive plotting off
plt.ioff()
################### input parameters #################################"
dir2 = '../results/aerosol-platt-nh3/aero/'

Nt= 5 * 60
#####################################################################


font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.pyplot.xticks(fontsize=14)

shape = [Nt]
poa1=np.zeros(shape)
poa2=np.zeros(shape)
poa3=np.zeros(shape)
soa1=np.zeros(shape)
soa2=np.zeros(shape)
soa3=np.zeros(shape)
soa4=np.zeros(shape)
soa5=np.zeros(shape)
soa6=np.zeros(shape)
inorg=np.zeros(shape)
isop=np.zeros(shape)

with open (dir2+'/PPOAlP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   poa1[j] = float(string_conc1[j])

with open (dir2+'/PPOAmP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   poa2[j] = float(string_conc1[j])

with open (dir2+'/PPOAhP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   poa3[j] = float(string_conc1[j])

with open (dir2+'/PSOAlP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa1[j] = float(string_conc1[j])

with open (dir2+'/PSOAmP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa2[j] = float(string_conc1[j])

with open (dir2+'/PSOAhP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa3[j] = float(string_conc1[j])

with open (dir2+'/PAnClP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa4[j] = float(string_conc1[j])

with open (dir2+'/PAnBlP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa5[j] = float(string_conc1[j])

with open (dir2+'/PAnBmP_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   soa6[j] = float(string_conc1[j])

with open (dir2+'/PNO3_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   inorg[j] = float(string_conc1[j])

with open (dir2+'/PNH4_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   inorg[j] = inorg[j] + float(string_conc1[j])

with open (dir2+'/PSO4_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   inorg[j] = inorg[j] + float(string_conc1[j])

with open (dir2+'/PBiPER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop[j] = isop[j] + float(string_conc1[j])

with open (dir2+'/PBiDER_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop[j] = isop[j] + float(string_conc1[j])

with open (dir2+'/PBiMT_1.txt') as finit :
	string_conc1 = finit.read().splitlines()
        for j in range(Nt) :
       	   isop[j] = isop[j] + float(string_conc1[j])

ppoa = poa1 + poa2 + poa3
totsoa = poa1 + poa2 + poa3 + soa1 + soa2 + soa3 + soa4 + soa5 + soa6 + inorg + isop
ssoa = soa1+soa2 + soa3
hsoa = soa5 + soa6

time=range(60,60*(Nt+1),60)
time=numpy.array(time)
time=time/3600.0

fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
fig = plt.figure(1,figsize = (15,15))
plt.ylim(0.,450.)

plt.plot(time, totsoa,'k-',label='total')
plt.plot(time, hsoa,color='navy',linestyle='dotted',label='SOA from VOC oxid. - high NOx')
plt.plot(time, soa4,'b-.',label='SOA from VOC oxid. - low NOx')
plt.plot(time, ssoa,'c--',label='SOA from SVOC oxid.')
plt.plot(time, ppoa,'g-',label='Primary SVOC')
plt.plot(time, inorg,'y-',label='Inorganics')
plt.ylabel("$\mu g.m^{-3}$", fontsize=18)
plt.xlabel("time [h]", fontsize=18)
plt.legend(loc=2)
plt.savefig("platt-particles-nh3.png",bbox_inches='tight')
plt.tight_layout()
plt.close(fig)
